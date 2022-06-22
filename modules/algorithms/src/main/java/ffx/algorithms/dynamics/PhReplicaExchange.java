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
import ffx.algorithms.Terminatable;
import ffx.potential.extended.ExtendedSystem;
import org.apache.commons.math3.util.FastMath;

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
public class PhReplicaExchange implements Terminatable {

  private static final Logger logger = Logger.getLogger(PhReplicaExchange.class.getName());
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
  /**
   * Each parameter array is wrapped inside a Parallel Java DoubleBuf for the All-Gather
   * communication calls.
   */
  private final DoubleBuf[] parametersBuf;

  private final MolecularDynamics replica;
  private boolean done = false;
  private boolean terminate = false;
  private final double[] myParameters;
  private final DoubleBuf myParametersBuf;

  private final int[] pH2Rank;
  private final int[] rank2Ph;
  private double[] pHScale;
  private final int[] pHAcceptedCount;
  private final int[] rankAcceptedCount;
  private final int[] pHTrialCount;
  private final double temp;

  private final ExtendedSystem extendedSystem;
  private final double pH;
  private double gapSize;

  /**
   * ReplicaExchange constructor.
   *
   * @param molecularDynamics a {@link MolecularDynamics} object.
   * @param pH pH = pKa <-- will be changed from this initial value
   * @param extendedSystem extended system attached to this process
   */
  public PhReplicaExchange(
      MolecularDynamics molecularDynamics, double pH, double temp, ExtendedSystem extendedSystem) {

    this.replica = molecularDynamics;
    this.temp = temp;
    this.extendedSystem = extendedSystem;
    this.pH = pH;

    // Set up the Replica Exchange communication variables for Parallel Java communication between
    // nodes.
    world = Comm.world();

    // Number of processes is equal to the number of replicas.
    int numProc = world.size();
    rank = world.rank();

    nReplicas = numProc;
    pHScale = new double[nReplicas];
    pH2Rank = new int[nReplicas];
    rank2Ph = new int[nReplicas];
    pHAcceptedCount = new int[nReplicas];
    rankAcceptedCount = new int [nReplicas];
    pHTrialCount = new int[nReplicas];

    setEvenSpacePhLadder(5, 12);

    random = new Random();
    random.setSeed(0);

    // Create arrays to store the parameters of all processes.
    parameters = new double[nReplicas][4]; //
    parametersBuf = new DoubleBuf[nReplicas];
    for (int i = 0; i < nReplicas; i++) {
      parametersBuf[i] = DoubleBuf.buffer(parameters[i]);
    }

    // A convenience reference to the parameters of this process are updated
    // during communication calls.
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
  public void sample(
      int cycles, long nSteps, double timeStep, double printInterval, double saveInterval) {
    done = false;
    terminate = false;
    for (int i = 0; i < cycles; i++) {
      // Check for termination request.
      if (terminate) {
        done = true;
        break;
      }
      dynamics(nSteps, timeStep, printInterval, saveInterval);
      logger.info(String.format(" Applying exchange condition for cycle %d.", i));
      exchange();

      logger.info(" Setting rank " + rank + " esv to pH " + pHScale[rank2Ph[rank]]);
    }

    logger.info("Replica Exchange Complete");
  }

  /**
   * setEvenPhLadder.
   *
   * @param lowpH a double.
   */
  public void setEvenSpacePhLadder(double lowpH, double highpH){
    gapSize = (highpH - lowpH) / nReplicas;
    for(int i = 0; i < nReplicas; i++){
      pHScale[i] = lowpH + i * gapSize;
      rank2Ph[i] = i;
      pH2Rank[i] = i;
    }
  }

/*
  public void setScaledPhLadder(double pHThatEqualsPKa){

    if(pHThatEqualsPKa > 2 && pHThatEqualsPKa < 12) {
      int lower = (int) (nReplicas * (pHThatEqualsPKa / 14.0)) - 1;
      int higher = (int) ( nReplicas * ((14 - pHThatEqualsPKa) / 14.0));

      for (int i = lower; i > 1; i--) {
        pHScale[lower - i] = pHThatEqualsPKa - ()
      }

      pHScale[lower] = pHThatEqualsPKa;

    }
    else {
      setEvenSpacePhLadder(pHThatEqualsPKa - 2);
    }
  }
 */

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
    for (int pH = 0; pH < nReplicas - 1; pH++) {

      // Ranks for pH A and B
      int rankA = pH2Rank[pH];
      int rankB = pH2Rank[pH + 1];

      // Load pH, beta and energy for each rank.
      double pHA = parameters[rankA][0];
      double pHB = parameters[rankB][0];
      double beta = KCAL_TO_GRAM_ANG2_PER_PS2 / (temp * kB);
      double acidostatA = parameters[rankA][2];
      double acidostatB = parameters[rankB][2];
      double acidostatAatB = parameters[rankA][3]; // acidostat of rankA evaluated at the pH of rankB
      double acidostatBatA = parameters[rankB][1];

      // Compute the change in energy over kT (E/kT) for the Metropolis criteria.
      logger.info("pHA = " + pHA + "" + acidostatA + "" + acidostatAatB);
      logger.info("pHB = " + pHB + "" + acidostatB + "" + acidostatBatA);
      logger.info("exp(" + beta + " * ((" + acidostatAatB + " - " + acidostatBatA + ") - (" + acidostatA + " + " + acidostatB + ")))");
      double deltaE = beta * ((acidostatAatB - acidostatBatA) - (acidostatA + acidostatB));

      //Count the number of trials for each temp
      pHTrialCount[pH]++;

      // If the Metropolis criteria is satisfied, do the switch.
      if (deltaE < 0.0 || random.nextDouble() < exp(-deltaE)) {
        pHAcceptedCount[pH]++;
        double pHAcceptance = pHAcceptedCount[pH] * 100.0 / (pHTrialCount[pH]);

        double rankAcceptance;
        if (pHA < pHB){
          rankAcceptedCount[rankA]++;
          rankAcceptance = rankAcceptedCount[rankA] * 100.0 / (pHTrialCount[pH]);
        } else {
          rankAcceptedCount[rankB]++;
          rankAcceptance = rankAcceptedCount[rankB] * 100.0 / (pHTrialCount[pH + 1]);
        }

        // Swap pH and energy values.
        parameters[rankA][0] = pHB;
        parameters[rankB][0] = pHA;

        // Map temperatures to process ranks.
        pH2Rank[pH] = rankB;
        pH2Rank[pH + 1] = rankA;

        // Map ranks to temperatures.
        rank2Ph[rankA] = pH + 1;
        rank2Ph[rankB] = pH;

        logger.info(
            String.format(
                " RepEx accepted. pH Accept: (%5.1f%%) Rank Accept: (%5.1f%%) for %6.2f (%d) and %6.2f (%d) for dE=%10.4f.",
                pHAcceptance, rankAcceptance, pHA, rankA, pHB, rankB, deltaE));
      } else {
        double tempAcceptance = pHAcceptedCount[pH] * 100.0 / (pHTrialCount[pH]);
        double rankAcceptance = rankAcceptedCount[pH] * 100.0 / (pHTrialCount[pH]);
        logger.info(
            String.format(
                " RepEx rejected (%5.1f%%) (f%5.1f%%) for %6.2f (%d) and %6.2f (%d) for dE=%10.4f.",
                tempAcceptance, rankAcceptance, pHA, rankA, pHB, rankB, deltaE));
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
  private void dynamics(
      final long nSteps,
      final double timeStep,
      final double printInterval,
      final double saveInterval) {

    int i = rank2Ph[rank];

    extendedSystem.setConstantPh(pHScale[i]);

    // Start this processes MolecularDynamics instance sampling.
    boolean initVelocities = true;
    replica.dynamic(
        nSteps, timeStep, printInterval, saveInterval, temp, initVelocities, null);

    // Update this ranks' parameter array to be consistent with the dynamics.
    myParameters[0] = pHScale[i];
    myParameters[2] = extendedSystem.getBiasEnergy();

    // Evaluate acidostat of ES at different pHs
    logger.info("Evaluating rank " + rank + " (originally at pH " + myParameters[0] + ") at pH " + (myParameters[0] - gapSize) + " and pH " + (myParameters[0] + gapSize));
    extendedSystem.setConstantPh(myParameters[0] - gapSize);
    myParameters[1] = extendedSystem.getBiasEnergy();

    extendedSystem.setConstantPh(myParameters[0] - gapSize);
    myParameters[3] = extendedSystem.getBiasEnergy();

    extendedSystem.setConstantPh(myParameters[0]);

    // Gather all parameters from the other processes.
    try {
      world.allGather(myParametersBuf, parametersBuf);
    } catch (IOException ex) {
      String message = " Replica Exchange allGather failed.";
      logger.log(Level.SEVERE, message, ex);
    }
  }
}
