//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.algorithms;

import edu.rit.mp.DoubleBuf;
import edu.rit.mp.IntegerBuf;
import edu.rit.pj.Comm;
import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parsers.SystemFilter;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.io.File;
import java.util.logging.Logger;

import static java.lang.String.format;
import static java.lang.System.arraycopy;

/**
 * The ParallelStateEnergy class evaluates the energy of a system at different lambda values.
 *
 * @since 1.0
 * @author Michael J. Schnieders
 */
public class ParallelStateEnergy {

  private static final Logger logger = Logger.getLogger(ParallelStateEnergy.class.getName());

  /**
   * Parallel Java world communicator.
   */
  private Comm world;
  /**
   * If false, do not use MPI communication.
   */
  private final boolean useMPI;
  /**
   * Number of processes.
   */
  private int numProc;
  /**
   * Rank of this process.
   */
  private int rank;
  /**
   * Number of states.
   */
  private int nStates;
  /**
   * Lambda value for each state.
   */
  private double[] lambdaValues;
  /**
   * The amount of work based on windows for each process.
   */
  private int statesPerProcess;

  /**
   * Number of samples for each state. The array is of size [statesPerProcess].
   */
  private final int[] nSamples;
  /**
   * The energy from evaluating at L - dL. The array is of size [statesPerProcess][snapshots].
   */
  private final double[][] energiesLowPJ;
  /**
   * The energy from evaluating at L. The array is of size [statesPerProcess][snapshots].
   */
  private final double[][] energiesAtPJ;
  /**
   * The energy from evaluating at L + dL. The array is of size [statesPerProcess][snapshots].
   */
  private final double[][] energiesHighPJ;
  /**
   * The volume of each snapshot. The array is of size [statesPerProcess][snapshots].
   */
  private final double[][] volumePJ;
  /**
   * Number of samples for each state. The array is of size [statesPerProcess].
   */
  private final IntegerBuf bufferNSamples;
  /**
   * Convenience reference for the DoubleBuf of this process.
   */
  private DoubleBuf bufferLow;
  /**
   * Convenience reference for the DoubleBuf of this process.
   */
  private DoubleBuf bufferAt;
  /**
   * Convenience reference for the DoubleBuf of this process.
   */
  private DoubleBuf bufferHigh;
  /**
   * Convenience reference for the DoubleBuf of this process.
   */
  private DoubleBuf bufferVolume;

  /**
   * The MolecularAssembly for each topology.
   */
  private final MolecularAssembly[] molecularAssemblies;
  /**
   * The SystemFilter for each topology.
   */
  private final SystemFilter[] openers;
  /**
   * The potential to evaluate.
   */
  private final Potential potential;
  /**
   * The full file paths for each state.
   */
  private final String[][] fullFilePaths;

  /**
   * The ParallelEnergy constructor.
   *
   * @param nStates             The number of states.
   * @param lambdaValues        The lambda values.
   * @param molecularAssemblies The molecular assemblies.
   * @param potential           The potential to evaluate.
   * @param fullFilePaths       The full file paths for each state.
   * @param openers             The system filters.
   */
  public ParallelStateEnergy(int nStates, double[] lambdaValues,
                             MolecularAssembly[] molecularAssemblies, Potential potential,
                             String[][] fullFilePaths, SystemFilter[] openers) {

    this.nStates = nStates;
    this.lambdaValues = lambdaValues;
    this.molecularAssemblies = molecularAssemblies;
    this.potential = potential;
    this.fullFilePaths = fullFilePaths;
    this.openers = openers;

    // Default to a single process that processes all states.
    world = Comm.world();
    numProc = 1;
    rank = 0;
    statesPerProcess = nStates;

    // Determine if use of PJ is specified.
    CompositeConfiguration properties = molecularAssemblies[0].getProperties();
    useMPI = properties.getBoolean("pj.use.mpi", true);
    if (useMPI) {
      // Number of processes.
      numProc = world.size();
      // Each processor gets its own rank.
      rank = world.rank();

      // Padding of the target array size (inner loop limit) is for parallelization.
      // Target states are parallelized over available nodes.
      // For example, if numProc = 8 and nStates = 12, then paddednWindows = 16.
      int extra = nStates % numProc;
      int paddednWindows = nStates;
      if (extra != 0) {
        paddednWindows = nStates - extra + numProc;
      }
      statesPerProcess = paddednWindows / numProc;

      if (numProc > 1) {
        logger.fine(format(" Number of MPI Processes:  %d", numProc));
        logger.fine(format(" Rank of this MPI Process: %d", rank));
        logger.fine(format(" States per process per row: %d", statesPerProcess));
      }
    }

    // Initialize arrays for storing energy values.
    nSamples = new int[statesPerProcess];
    energiesLowPJ = new double[statesPerProcess][];
    energiesAtPJ = new double[statesPerProcess][];
    energiesHighPJ = new double[statesPerProcess][];
    volumePJ = new double[statesPerProcess][];
    bufferNSamples = IntegerBuf.buffer(nSamples);
    bufferLow = DoubleBuf.buffer(energiesLowPJ);
    bufferAt = DoubleBuf.buffer(energiesAtPJ);
    bufferHigh = DoubleBuf.buffer(energiesHighPJ);
    bufferVolume = DoubleBuf.buffer(volumePJ);
  }

  /**
   * Get the rank of this process.
   *
   * @return The rank.
   */
  public int getRank() {
    return rank;
  }

  /**
   * Evaluate the energies for each state.
   *
   * @param energiesLow  The energy from evaluating at L - dL.
   * @param energiesAt   The energy from evaluating at L.
   * @param energiesHigh The energy from evaluating at L + dL.
   * @param volume       The volume of each snapshot.
   */
  public void evaluateStates(double[][] energiesLow,
                             double[][] energiesAt,
                             double[][] energiesHigh,
                             double[][] volume) {
    double[] currentLambdas;
    double[][] energy;
    int nCurrLambdas;
    for (int state = 0; state < nStates; state++) {
      if (state == 0) {
        currentLambdas = new double[2];
        currentLambdas[0] = lambdaValues[state];
        currentLambdas[1] = lambdaValues[state + 1];
      } else if (state == nStates - 1) {
        currentLambdas = new double[2];
        currentLambdas[0] = lambdaValues[state - 1];
        currentLambdas[1] = lambdaValues[state];
      } else {
        currentLambdas = new double[3];
        currentLambdas[0] = lambdaValues[state - 1];
        currentLambdas[1] = lambdaValues[state];
        currentLambdas[2] = lambdaValues[state + 1];
      }
      nCurrLambdas = currentLambdas.length;
      energy = new double[nCurrLambdas][];
      evaluateEnergies(state, currentLambdas, energy, fullFilePaths);
    }

    gatherAllValues(energiesLow, energiesAt, energiesHigh, volume);
  }


  /**
   * Evaluate the energies for a given state.
   *
   * @param state         The state.
   * @param lambdaValues  The current lambda values.
   * @param energy        The energy values.
   * @param fullFilePaths The full file paths.
   */
  private void evaluateEnergies(int state, double[] lambdaValues,
                                double[][] energy, String[][] fullFilePaths) {
    if (state % numProc == rank) {
      int workItem = state / numProc;
      double[] vol = getEnergyForLambdas(lambdaValues, fullFilePaths[state], energy);
      int len = energy[0].length;
      nSamples[workItem] = len;
      energiesLowPJ[workItem] = new double[len];
      energiesAtPJ[workItem] = new double[len];
      energiesHighPJ[workItem] = new double[len];
      volumePJ[workItem] = new double[len];
      if (state == 0) {
        arraycopy(energy[0], 0, energiesAtPJ[workItem], 0, len);
        arraycopy(energy[1], 0, energiesHighPJ[workItem], 0, len);
      } else if (state < nStates - 1) {
        arraycopy(energy[0], 0, energiesLowPJ[workItem], 0, len);
        arraycopy(energy[1], 0, energiesAtPJ[workItem], 0, len);
        arraycopy(energy[2], 0, energiesHighPJ[workItem], 0, len);
      } else if (state == nStates - 1) {
        arraycopy(energy[0], 0, energiesLowPJ[workItem], 0, len);
        arraycopy(energy[1], 0, energiesAtPJ[workItem], 0, len);
      }
      if (vol != null) {
        arraycopy(vol, 0, volumePJ[workItem], 0, len);
      }
    }
  }

  /**
   * Gather all energy values from all nodes.
   * This method calls <code>world.gather</code> to collect numProc values.
   *
   * @param energiesLow  The energy from evaluating at L - dL.
   * @param energiesAt   The energy from evaluating at L.
   * @param energiesHigh The energy from evaluating at L + dL.
   * @param volume       The volume of each snapshot.
   */
  private void gatherAllValues(double[][] energiesLow,
                               double[][] energiesAt,
                               double[][] energiesHigh,
                               double[][] volume) {
    if (useMPI) {
      try {
        if (rank != 0) {
          // Send all results to node 0.
          world.send(0, bufferNSamples);
          for (int workItem = 0; workItem < statesPerProcess; workItem++) {
            bufferLow = DoubleBuf.buffer(energiesLowPJ[workItem]);
            bufferAt = DoubleBuf.buffer(energiesAtPJ[workItem]);
            bufferHigh = DoubleBuf.buffer(energiesHighPJ[workItem]);
            bufferVolume = DoubleBuf.buffer(volumePJ[workItem]);
            world.send(0, bufferLow);
            world.send(0, bufferAt);
            world.send(0, bufferHigh);
            world.send(0, bufferVolume);
          }
        } else {
          for (int proc = 0; proc < numProc; proc++) {
            // Receive all results from another node.
            if (proc > 0) {
              // Receive the number of samples for each state.
              world.receive(proc, bufferNSamples);
            }
            // Store results in the appropriate arrays.
            for (int workItem = 0; workItem < statesPerProcess; workItem++) {
              final int state = numProc * workItem + proc;
              final int nSnapshots = nSamples[workItem];
              // Do not include padded results.
              if (state < nStates) {
                if (proc > 0) {
                  // Ensure array size.
                  updateMemory(workItem, nSnapshots);
                  world.receive(proc, bufferLow);
                  world.receive(proc, bufferAt);
                  world.receive(proc, bufferHigh);
                  world.receive(proc, bufferVolume);
                }
                energiesLow[state] = new double[nSnapshots];
                energiesAt[state] = new double[nSnapshots];
                energiesHigh[state] = new double[nSnapshots];
                volume[state] = new double[nSnapshots];
                arraycopy(energiesLowPJ[workItem], 0, energiesLow[state], 0, nSnapshots);
                arraycopy(energiesAtPJ[workItem], 0, energiesAt[state], 0, nSnapshots);
                arraycopy(energiesHighPJ[workItem], 0, energiesHigh[state], 0, nSnapshots);
                arraycopy(volumePJ[workItem], 0, volume[state], 0, nSnapshots);
              }
            }
          }
        }
        // Wait for all processes to finish.
        world.barrier();
      } catch (Exception ex) {
        logger.severe(" Exception collecting energy values." + ex + Utilities.stackTraceToString(ex));
      }
    } else {
      for (int i = 0; i < nStates; i++) {
        int len = energiesAtPJ[rank].length;
        energiesLow[i] = new double[len];
        energiesAt[i] = new double[len];
        energiesHigh[i] = new double[len];
        volume[i] = new double[len];
        arraycopy(energiesLowPJ[i], 0, energiesLow[i], 0, len);
        arraycopy(energiesAtPJ[i], 0, energiesAt[i], 0, len);
        arraycopy(energiesHighPJ[i], 0, energiesHigh[i], 0, len);
        arraycopy(volumePJ[i], 0, volume[i], 0, len);
      }
    }
  }

  /**
   * Compute energy values for each lambda value.
   *
   * @param lambdaValues The lambda values.
   * @param arcFileName  The archive file names.
   * @param energy       The energy values.
   * @return The volume of each snapshot.
   */
  private double[] getEnergyForLambdas(double[] lambdaValues,
                                       String[] arcFileName, double[][] energy) {

    int numTopologies = molecularAssemblies.length;

    // Initialize the potential to use the correct archive files.
    StringBuilder sb = new StringBuilder("\n");
    for (int j = 0; j < numTopologies; j++) {
      File archiveFile = new File(arcFileName[j]);
      openers[j].setFile(archiveFile);
      molecularAssemblies[j].setFile(archiveFile);
      sb.append(format(" Evaluating energies for file: %s\n", arcFileName[j]));
    }
    sb.append("\n");
    logger.info(sb.toString());

    int nSnapshots = openers[0].countNumModels();
    double[] x = new double[potential.getNumberOfVariables()];
    double[] vol = new double[nSnapshots];
    int nLambdas = lambdaValues.length;
    for (int k = 0; k < nLambdas; k++) {
      energy[k] = new double[nSnapshots];
    }

    LambdaInterface linter1 = (LambdaInterface) potential;

    int endWindow = nStates - 1;
    String endWindows = endWindow + File.separator;

    if (arcFileName[0].contains(endWindows)) {
      logger.info(format(" %s     %s   %s     %s   %s ", "Snapshot", "Lambda Low",
          "Energy Low", "Lambda At", "Energy At"));
    } else if (arcFileName[0].contains("0/")) {
      logger.info(format(" %s     %s   %s     %s   %s ", "Snapshot", "Lambda At",
          "Energy At", "Lambda High", "Energy High"));
    } else {
      logger.info(format(" %s     %s   %s     %s   %s     %s   %s ", "Snapshot", "Lambda Low",
          "Energy Low", "Lambda At", "Energy At", "Lambda High", "Energy High"));
    }

    for (int i = 0; i < nSnapshots; i++) {
      boolean resetPosition = (i == 0);
      int nOpeners = openers.length;
      for (int n = 0; n < nOpeners; n++) {
        openers[n].readNext(resetPosition, false);
      }

      x = potential.getCoordinates(x);
      nLambdas = lambdaValues.length;
      for (int k = 0; k < nLambdas; k++) {
        double lambda = lambdaValues[k];
        linter1.setLambda(lambda);
        energy[k][i] = potential.energy(x, false);
      }

      if (nLambdas == 2) {
        logger.info(format(" %8d     %6.3f   %14.4f     %6.3f   %14.4f ", i + 1,
            lambdaValues[0], energy[0][i], lambdaValues[1], energy[1][i]));
      } else {
        logger.info(format(" %8d     %6.3f   %14.4f     %6.3f   %14.4f     %6.3f   %14.4f ", i + 1,
            lambdaValues[0], energy[0][i], lambdaValues[1], energy[1][i], lambdaValues[2], energy[2][i]));
      }

      Crystal unitCell;
      if (potential instanceof CrystalPotential) {
        unitCell = ((CrystalPotential) potential).getCrystal().getUnitCell();
      } else {
        unitCell = molecularAssemblies[0].getCrystal().getUnitCell();
      }

      if (!unitCell.aperiodic()) {
        int nSymm = unitCell.getNumSymOps();
        vol[i] = unitCell.volume / nSymm;
      }
    }

    return vol;
  }

  /**
   * Update the memory to receive energy values given the number of snapshots.
   *
   * @param workItem   The work item to update.
   * @param nSnapshots The number of snapshots.
   */
  private void updateMemory(int workItem, int nSnapshots) {
    if (energiesAtPJ[workItem].length < nSnapshots) {
      energiesLowPJ[workItem] = new double[nSnapshots];
      energiesAtPJ[workItem] = new double[nSnapshots];
      energiesHighPJ[workItem] = new double[nSnapshots];
      volumePJ[workItem] = new double[nSnapshots];
    }
    bufferLow = DoubleBuf.buffer(energiesLowPJ[workItem]);
    bufferAt = DoubleBuf.buffer(energiesAtPJ[workItem]);
    bufferHigh = DoubleBuf.buffer(energiesHighPJ[workItem]);
    bufferVolume = DoubleBuf.buffer(volumePJ[workItem]);
  }
}
