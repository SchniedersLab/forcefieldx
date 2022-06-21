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
package ffx.algorithms.dynamics;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import ffx.potential.extended.ExtendedSystem;

import java.io.IOException;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.kB;
import static org.apache.commons.math3.util.FastMath.exp;

/**
 * The ReplicaExchange implements temperature and lambda replica exchange methods.
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 * @since 1.0
 */
public class TempPhReplicaExchange implements Terminatable {

  private static final Logger logger = Logger.getLogger(TempPhReplicaExchange.class.getName());
  private final int nReplicas;
  private final Random random;
  /** Parallel Java world communicator. */
  private final Comm world;
  /** Rank of this process. */
  private final int rank;
  /**
   * The parameters array stores communicated parameters for each process
   * (i.e. each RepEx system).
   * Currently the array is of size [number of Processes][2].
   */
  private final double[][] parameters;
  private final double[] myParameters;
  /**
   * Each parameter array is wrapped inside a Parallel Java DoubleBuf for the All-Gather
   * communication calls.
   */
  private final DoubleBuf[] parametersBuf;
  private final DoubleBuf myParametersBuf;
  private final MolecularDynamics replica;
  private boolean done = false;
  private boolean terminate = false;
  private final AlgorithmListener listener;
  private final ExtendedSystem extendedSystem;
  private final double pH, temp;
  private final int pHRows, tempColumns;
  private final int[][] rankArray;
  private final double[] tempScale, pHScale;
  private final int[] pHTrialCount, pHAcceptedCount, tempTrialCount, tempAcceptedCount, rankAcceptedCount;


  /**
   * ReplicaExchange constructor.
   *
   * @param molecularDynamics a {@link MolecularDynamics} object.
   * @param listener a {@link AlgorithmListener} object.
   * @param pH pH = pKa <-- will be changed from this initial value
   * @param extendedSystem extended system attached to this process
   */
  public TempPhReplicaExchange(
      MolecularDynamics molecularDynamics, AlgorithmListener listener, double pH, int pHRows , double temp, int tempColumns, double exponent, ExtendedSystem extendedSystem) {

    // Set up the Replica Exchange communication variables for Parallel Java communication between nodes.
    this.world = Comm.world();

    // Number of processes is equal to the number of replicas.
    this.nReplicas = world.size();
    this.rank = world.rank();
    this.rankArray = new int[pHRows][tempColumns];

    if(nReplicas != pHRows * tempColumns) {
      logger.severe("nReplicas does not match pHRows * tempColumns. Program will go out of bounds.");
    }

    for (int i = 0; i < pHRows; i++) {
      for (int j = 0; j < tempColumns; j++){
        rankArray[i][j] = i*j + j;
      }
    }

    //Identical Seed for each process
    this.random = new Random();
    this.random.setSeed(0);

    // Initialized option variables
    this.replica = molecularDynamics;
    this.listener = listener;
    this.extendedSystem = extendedSystem;
    this.pH = pH;
    this.temp = temp;
    this.pHRows = pHRows;
    this.tempColumns = tempColumns;

    // Automatic variables based on options
    pHScale = new double[pHRows];
    tempScale = new double[tempColumns];

    pHTrialCount = new int[pHRows];
    pHAcceptedCount = new int[pHRows];
    tempTrialCount = new int[tempColumns];
    tempAcceptedCount = new int[tempColumns];
    rankAcceptedCount = new int[nReplicas];

    setEvenSpacePhLadder(2, 12);
    setExponentialTemperatureLadder(273, exponent);

    // Create arrays to store the parameters of all processes.
    parameters = new double[nReplicas][3]; // parameters[rank][index] | index = 0-pH 1-temp 2-energy
    parametersBuf = new DoubleBuf[nReplicas];
    for (int i = 0; i < nReplicas; i++) {
      parametersBuf[i] = DoubleBuf.buffer(parameters[i]);
    }

    // A convenience reference to the parameters of this process are updated during communication calls.
    myParameters = parameters[rank];
    myParametersBuf = parametersBuf[rank];
  }

  /**
   * Sample.
   *
   * @param cycles a int.
   * @param nSteps a int.
   * @param timeStep a double.
   * @param printInterval a double.
   * @param saveInterval a double.
   */
  public void sample(int cycles, long nSteps, double timeStep, double printInterval, double saveInterval) {
    done = false;
    terminate = false;
    for (int i = 0; i < cycles; i++) {
      // Check for termination request.
      if (terminate) {
        done = true;
        break;
      }
      dynamic(nSteps, timeStep, printInterval, saveInterval);
      logger.info(String.format(" Applying exchange condition for cycle %d.", i));
      exchange();

      //TODO: change pH of system before next dynamic cycle and after evaluating exchange();
      extendedSystem.setConstantPh(parameters[rank][0]);
      extendedSystem.writeLambdaHistogram();
    }
  }

  /**
   * setExponentialTemperatureLadder.
   *
   * @param lowTemperature a double.
   * @param exponent a double.
   */
  public void setExponentialTemperatureLadder(double lowTemperature, double exponent) {
    for (int i = 0; i < tempColumns; i++) {
      tempScale[i] = lowTemperature * exp(exponent * i);

      for(int j = 0; j < pHRows; j++){
        int index = i*j + j;
        parameters[index][1] = tempScale[i];
      }
    }
  }

  /**
   * setEvenPhLadder.
   *
   * @param lowpH a double.
   */
  public void setEvenSpacePhLadder(double lowpH, double highpH){
    for(int i = 0; i < pHRows; i++){
      pHScale[i] = lowpH + i * (int) ((highpH - lowpH ) / pHRows);

      for (int j = 0; j < tempColumns; j++){
        int index = i*j + i;
        parameters[index][0] = pHScale[i];
      }
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>This should be implemented as a blocking interrupt; when the method returns the <code>
   * Terminatable</code> algorithm has reached a clean termination point. For example, between
   * minimize or molecular dynamics steps.
   */
  @Override
  public void terminate() {
    terminate = true;
    while (!done) {
      synchronized (this) {
        try {
          wait(1);
        } catch (InterruptedException e) {
          logger.log(Level.WARNING, "Exception terminating replica exchange.\n", e);
        }
      }
    }
  }

  /** All processes complete the exchanges identically given the same Random number seed. */
  private void exchange() {
    // 2 M.C. trials per pH (except for those at the ends of the ladder).
    // Loop over pH scale
    for (int pH = 0; pH < pHScale.length - 1; pH++) {
      for (int temperature = 0; temperature < tempScale.length - 1; temperature++) {
        for(int i = 0; i < 3; i++) {
          // Ranks for index A and B
          int rankA = rankArray[pH][temperature];
          int rankB = -1;

          switch (i) {
            case 2:
              rankB = rankArray[pH][temperature + 1];
              break;
            case 1:
              rankB = rankArray[pH + 1][temperature];
              break;
            case 0:
              rankB = rankArray[pH + 1][temperature + 1];
              break;
          }

          // Load pH, temp, beta, and energy for each rank.
          double pHA = parameters[rankA][0];
          double pHB = parameters[rankB][0];

          double tempA = parameters[rankA][1];
          double tempB = parameters[rankB][1];

          double energyA = parameters[rankA][0];
          double energyB = parameters[rankB][0];

          double betaA = KCAL_TO_GRAM_ANG2_PER_PS2 / (tempA * kB);
          double betaB = KCAL_TO_GRAM_ANG2_PER_PS2 / (tempB * kB);


          // Compute the change in energy over kT (E/kT) for the Metropolis criteria.
          double deltaE = (energyA - energyB) * (betaB - betaA);

          //Count the number of trials for each temp
          switch (i) {
            case 2:
              tempTrialCount[temperature]++;
              break;
            case 1:
              pHTrialCount[pH]++;
              break;
            case 0:
              tempTrialCount[temperature]++;
              pHTrialCount[pH]++;
              break;
          }

          // If the Metropolis criteria is satisfied, do the switch.
          if (deltaE < 0.0 || random.nextDouble() < exp(-deltaE)) {

            // Swap Ranks
            rankArray[pH][temperature] = rankB;
            switch (i) {
              case 2:
                tempAcceptedCount[temperature]++;
                rankArray[pH][temperature + 1] = rankA;
                break;
              case 1:
                pHAcceptedCount[pH]++;
                rankArray[pH + 1][temperature] = rankA;
                break;
              case 0:
                rankArray[pH + 1][temperature + 1] = rankA;
                tempAcceptedCount[temperature]++;
                pHAcceptedCount[pH]++;
                break;
            }


            double tempAcceptance = tempAcceptedCount[pH] * 100.0 / (tempTrialCount[pH]);
            double pHAcceptance = pHAcceptedCount[pH] * 100.0 / (pHTrialCount[pH]);

            double rankAcceptance = -1;

            switch (i) {
              case 2:
                if (pHA < pHB) {
                  rankAcceptedCount[rankA]++;
                  rankAcceptance = rankAcceptedCount[rankA] * 100.0 / (pHTrialCount[pH]);
                } else {
                  rankAcceptedCount[rankB]++;
                  rankAcceptance = rankAcceptedCount[rankB] * 100.0 / (pHTrialCount[pH + 1]);
                }
                break;
              case 1:
                if (tempA < tempB) {
                  rankAcceptedCount[rankA]++;
                  rankAcceptance = rankAcceptedCount[rankA] * 100.0 / (tempTrialCount[temperature]);
                } else {
                  rankAcceptedCount[rankB]++;
                  rankAcceptance = rankAcceptedCount[rankB] * 100.0 / (tempTrialCount[temperature + 1]);
                }
                break;
              case 0:

                if (pHA < pHB) {
                  rankAcceptedCount[rankA]++;
                } else {
                  rankAcceptedCount[rankB]++;
                }

                if (tempA < tempB) {
                  rankAcceptedCount[rankA]++;
                  rankAcceptance = rankAcceptedCount[rankA] * 100.0 / (tempTrialCount[temperature]);
                } else {
                  rankAcceptedCount[rankB]++;
                  rankAcceptance = rankAcceptedCount[rankB] * 100.0 / (tempTrialCount[temperature + 1]);
                }
                break;
            }

            // Swap pH and energy and temperature values.
            parameters[rankA][0] = pHB;
            parameters[rankB][0] = pHA;

            parameters[rankA][1] = energyB;
            parameters[rankB][1] = energyA;

            parameters[rankA][2] = tempB;
            parameters[rankB][2] = tempA;

            switch (i) {
              case 0:
                logger.info("Accepted: rank[i][j] --> rank[i][j + 1]");
                break;
              case 1:
                logger.info("Accepted: rank[i][j] --> rank[i + 1][j]");
                break;
              case 2:
                logger.info("Accepted: rank[i][j] --> rank[i + 1][j + 1]");
                break;
            }

            logger.info(
                    String.format(
                            " RepEx accepted (%5.1f%%) (%5.1f%%) for %6.2f (%d) and %6.2f (%d) for dE=%10.4f.",
                            tempAcceptance, rankAcceptance, pHA, rankA, pHB, rankB, deltaE));

            logger.info(
                    String.format(
                            " RepEx accepted (%5.1f%%) (%5.1f%%) for %6.2f (%d) and %6.2f (%d) for dE=%10.4f.",
                            pHAcceptance, rankAcceptance, pHA, rankA, pHB, rankB, deltaE));

            break;

          } else {
            double tempAcceptance = pHAcceptedCount[pH] * 100.0 / (pHTrialCount[pH]);
            double rankAcceptance = rankAcceptedCount[pH] * 100.0 / (pHTrialCount[pH]);

            switch (i) {
              case 0:
                logger.info("Rejected: rank[i][j] --> rank[i][j + 1]");
                break;
              case 1:
                logger.info("Rejected: rank[i][j] --> rank[i + 1][j]");
                break;
              case 2:
                logger.info("Rejected: rank[i][j] --> rank[i + 1][j + 1]");
                break;
            }

            logger.info(
                    String.format(
                            " RepEx rejected (%5.1f%%) (f%5.1f%%) for %6.2f (%d) and %6.2f (%d) for dE=%10.4f.",
                            tempAcceptance, rankAcceptance, pHA, rankA, pHB, rankB, deltaE));

          }
        }
      }
    }
  }

  /**
   * Blocking dynamic steps: when this method returns each replica has completed the requested
   * number of steps.
   *
   * @param nSteps the number of time steps.
   * @param timeStep the time step.
   * @param printInterval the number of steps between loggging updates.
   * @param saveInterval the number of steps between saving snapshots.
   */
  private void dynamic(
      final long nSteps,
      final double timeStep,
      final double printInterval,
      final double saveInterval) {

    // Start this processes MolecularDynamics instance sampling.
    boolean initVelocities = true;
    replica.dynamic(
        nSteps, timeStep, printInterval, saveInterval, myParameters[2], initVelocities, null);


    // Update this ranks' parameter array to be consistent with the dynamics.
    myParameters[1] = replica.currentPotentialEnergy;
    for (int i = 0; i < pHRows; i++){
      for (int j = 0; j < tempColumns; j++){
        if(rankArray[i][j] == rank){
          myParameters[0] = pHScale[i];
          myParameters[2] = tempScale[j];
        }
      }
    }

    // Gather all parameters from the other processes.
    try {
      world.allGather(myParametersBuf, parametersBuf);
    } catch (IOException ex) {
      String message = " Replica Exchange allGather failed.";
      logger.log(Level.SEVERE, message, ex);
    }
  }
}
