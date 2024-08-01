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
package ffx.algorithms.dynamics;

import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.kB;
import static org.apache.commons.math3.util.FastMath.exp;

import edu.rit.mp.DoubleBuf;
import edu.rit.mp.IntegerBuf;
import edu.rit.mp.buf.IntegerMatrixBuf_1;
import edu.rit.pj.Comm;
import ffx.algorithms.Terminatable;
import ffx.numerics.Potential;
import ffx.potential.extended.ExtendedSystem;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.potential.parameters.TitrationUtils;
import jline.internal.Nullable;
import org.apache.commons.io.FilenameUtils;

/**
 * The PhReplicaExchange implements pH replica exchange. Adapted from "ReplicaExchange.java" by
 * Timothy D. Fenn and Michael J. Schnieders
 *
 * @author Matthew Speranza and Andrew Thiel
 * @since 1.0
 */
public class PhReplicaExchange implements Terminatable {

  private static final Logger logger = Logger.getLogger(PhReplicaExchange.class.getName());
  private final MolecularDynamics replica;
  private final MolecularDynamicsOpenMM openMM;
  private final Potential potential;
  private final ExtendedSystem extendedSystem;
  private final Random random;
  /** Parallel Java world communicator. */
  private final Comm world;
  private File dyn, esv, esvBackup, dynBackup;
  private ArrayList<Double> readPhScale;
  /**
   * Each parameter array is wrapped inside a Parallel Java DoubleBuf for the All-Gather
   * communication calls.
   */
  private final DoubleBuf[] parametersBuf;
  private final DoubleBuf myParametersBuf;
  /**
   * The parameters array stores communicated parameters for each process (i.e. each RepEx system).
   * Currently, the array is of size [number of Processes][2].
   */
  private final double[][] parameters;
  private final int[][][] parametersHis;
  private final IntegerBuf[] parametersHisBuf;
  private final int[] pH2Rank, rank2Ph, pHAcceptedCount, pHTrialCount;
  /** Rank of this process. */
  private final int rank, nReplicas;
  private boolean done = false, restart, terminate = false, backupNeeded = false;
  private final double[] myParameters, pHScale;
  private double[] x;
  private final double pH, temp;
  private int restartStep;

  /**
   * pHReplicaExchange constructor.
   *
   * @param molecularDynamics a {@link MolecularDynamics} object.
   * @param pH pH = pKa will be changed from this initial value
   * @param extendedSystem extended system attached to this process
   * @param pHLadder range of pH's that replicas are created for
   * @param temp temperature of replica
   */
  public PhReplicaExchange(MolecularDynamics molecularDynamics, File structureFile, double pH,
                           double[] pHLadder, double temp, ExtendedSystem extendedSystem, double[] x) {
    this(molecularDynamics, structureFile, pH, pHLadder, temp, extendedSystem, x, null, null);
  }

  /**
   * OpenMM cycled pHReplicaExchange constructor.
   *
   * @param molecularDynamics a {@link MolecularDynamics} object.
   * @param pH pH = pKa will be changed from this initial value
   * @param extendedSystem extended system attached to this process
   * @param pHLadder range of pH's that replicas are created for
   * @param temp temperature of replica
   * @param x array of coordinates
   * @param molecularDynamicsOpenMM for running config steps on GPU
   */
  public PhReplicaExchange(MolecularDynamics molecularDynamics, File structureFile, double pH,
                           double[] pHLadder, double temp, ExtendedSystem extendedSystem, @Nullable double[] x,
                           @Nullable MolecularDynamicsOpenMM molecularDynamicsOpenMM, @Nullable Potential potential) {

    this.replica = molecularDynamics;
    this.temp = temp;
    this.extendedSystem = extendedSystem;
    this.pH = pH;
    this.pHScale = pHLadder;
    this.x = x;
    this.openMM = molecularDynamicsOpenMM;
    this.potential = potential;

    // Set up the Replica Exchange communication variables for Parallel Java communication between nodes.
    world = Comm.world();

    // Number of processes is equal to the number of replicas.
    int numProc = world.size();
    rank = world.rank();

    nReplicas = numProc;
    pH2Rank = new int[nReplicas];
    rank2Ph = new int[nReplicas];
    pHAcceptedCount = new int[nReplicas];
    pHTrialCount = new int[nReplicas];

    // Init map
    for (int i = 0; i < nReplicas; i++) {
      rank2Ph[i] = i;
      pH2Rank[i] = i;
    }
    extendedSystem.setConstantPh(pHScale[rank]);

    random = new Random();
    random.setSeed(0);

    // Create arrays to store the parameters of all processes.
    parameters = new double[nReplicas][pHLadder.length + 1];
    parametersHis = new int[nReplicas][extendedSystem.getTitratingResidueList().size()][100];
    parametersBuf = new DoubleBuf[nReplicas];
    parametersHisBuf = new IntegerMatrixBuf_1[nReplicas];
    for (int i = 0; i < nReplicas; i++) {
      parametersBuf[i] = DoubleBuf.buffer(parameters[i]);
      parametersHisBuf[i] = IntegerMatrixBuf_1.buffer(parametersHis[i]);
    }

    // A convenience reference to the parameters of this process are updated during communication calls.
    myParameters = parameters[rank];
    myParametersBuf = parametersBuf[rank];

    // Set the ESV that this rank will share has correct values after a restart
    if (manageFilesAndRestart(structureFile)) {
      logger.info(" Startup Successful! ");
      extendedSystem.getESVHistogram(parametersHis[rank]);
    } else {
      logger.severe(" Unable to restart. Fix issues with restart files.");
    }
  }

  /**
   * This deals with anything having to do with restarts and backup files. Is identical between CPU
   * and GPU.
   *
   * @param structureFile pdb file
   * @return whether the restart was successful or not
   */
  private boolean manageFilesAndRestart(File structureFile) {
    // My directory info
    File parent = structureFile.getParentFile();
    File rankDir = new File(parent + File.separator + rank);
    restart = !rankDir.mkdir();

    // My files
    esv = new File(
            rankDir.getPath() + File.separator + FilenameUtils.removeExtension(structureFile.getName())
                    + ".esv");
    dyn = new File(
            rankDir.getPath() + File.separator + FilenameUtils.removeExtension((structureFile.getName()))
                    + ".dyn");
    esvBackup = new File(FilenameUtils.removeExtension(esv.getAbsolutePath()) + "_backup.esv");
    dynBackup = new File(FilenameUtils.removeExtension(dyn.getAbsolutePath()) + "_backup.dyn");

    // Set restart files - does not do any reading yet
    extendedSystem.setRestartFile(esv);
    replica.setFallbackDynFile(dyn);
    replica.setArchiveFiles(new File[] {new File(
            rankDir.getPath() + File.separator + FilenameUtils.removeExtension((structureFile.getName()))
                    + ".arc")});
    if (openMM != null) {
      openMM.setAutomaticWriteouts(false);
    }

    // Build restart
    restartStep = 0;
    if (restart) {
      logger.info(" Attempting to restart. ");
      readPhScale = new ArrayList<>();

      // Read normal restarts
      backupNeeded = false;
      if (checkForRestartFiles(parent, esv.getName()) && checkForRestartFiles(parent,
              dyn.getName())) {
        for (int i = 0; i < nReplicas; i++) {
          File checkThisESV = new File(
                  parent.getAbsolutePath() + File.separator + i + File.separator + esv.getName());
          backupNeeded = !readESV(checkThisESV, i);
          if (backupNeeded) {
            logger.warning(" Searching backup esv files.");
            break;
          }
        }
        if (!backupNeeded) {
          logger.info(" Reading in from: " + esv);
          extendedSystem.readESVInfoFrom(esv);
        }
      } else {
        logger.info(" Directories do not contain all of the correct restart files.");
        backupNeeded = true;
      }

      // Read backup restarts
      if (backupNeeded && checkForRestartFiles(parent, esvBackup.getName()) && checkForRestartFiles(
              parent, dynBackup.getName())) {
        for (int i = 0; i < nReplicas; i++) {
          File checkThisESV = new File(
                  parent.getAbsolutePath() + File.separator + i + File.separator + esvBackup.getName());
          if (!readESV(checkThisESV, i)) {
            logger.info(" Restart files & Backups are messed up.");
            return false;
          }
        }
        logger.info(" Reading in from: " + esvBackup);
        extendedSystem.readESVInfoFrom(esvBackup);
      } else if (backupNeeded) {
        logger.info(" Directories do not contain all of the correct backup restart files.");
        logger.warning(" No backup files found, treating as a fresh start.");
        restart = false;
        return true;
      }
    } else {
      logger.info(" Not a restart.");
      return true;
    }

    logger.info(" Rank " + rank + " staring at pH " + pHScale[rank2Ph[rank]]);
    extendedSystem.setConstantPh(pHScale[rank2Ph[rank]]);
    logger.info(" Rank " + rank + " starting hist:");
    extendedSystem.writeLambdaHistogram(true);

    return true;
  }

  /**
   * Reads an esv file to find the number of steps logged and if the pH matches the others in
   * readPhScale. (N.B. - This only checks the first esv in the system for uneven sums, since it is
   * assumed that all counts will be messed up if the first one is.)
   *
   * @param file to read
   * @return if the esv file of this name at this rank works with the others, or if i=0, then it sets
   *     the standards
   */
  private boolean readESV(File file, int i) {
    try (BufferedReader br = new BufferedReader(new FileReader(file))) {
      String data = br.readLine();
      while (data != null) {
        List<String> tokens = Arrays.asList(data.split(" +"));
        if (tokens.contains("pH:")) {
          double pHOfRankI = Double.parseDouble(tokens.get(tokens.indexOf("pH:") + 1));

          // Set up map
          for (int j = 0; j < pHScale.length; j++) {
            if (Math.abs(pHScale[j] / pHOfRankI - 1) < 0.00001) {
              rank2Ph[i] = j;
              pH2Rank[j] = i;
            }
          }

          // Check pH list for duplicates
          if (readPhScale.stream().anyMatch(d -> (Math.abs(d / pHOfRankI - 1) < 0.00001))) {
            logger.warning(" Duplicate pH value found in file: " + file.getAbsolutePath());
            readPhScale.clear();
            restartStep = 0;
            return false;
          } else {
            readPhScale.add(pHOfRankI);
          }

          // Find restart steps
          br.readLine();
          br.readLine();
          int sum = 0;
          for (int j = 0; j < 10; j++) { // 10 titr windows
            data = br.readLine().trim();
            tokens = List.of(data.split(" +"));
            for (int k = 0; k < 10; k++) { // 10 tautomer windows
              sum += Integer.parseInt(tokens.get(k + 1));
            }
          }

          if (restartStep == 0) {
            restartStep = sum;
          } else if (restartStep != sum) {
            logger.warning(" Restart received uneven sums from esv file: " + file.getAbsolutePath());
            logger.info(" Restart Step: " + restartStep + "    Found: " + sum);
            restartStep = 0;
            readPhScale.clear();
            return false;
          }
          break;
        }
        data = br.readLine();
      }
      return true;
    } catch (IOException | IndexOutOfBoundsException e) {
      logger.warning("Failed to read file: " + file.getAbsolutePath());
      return false;
    }
  }

  /**
   * Looks through directories to see if they all contain the right file
   *
   * @return whether they all contain the correct file
   */
  private boolean checkForRestartFiles(File parent, String fileNameWithExtension) {
    File searchFile;
    for (int i = 0; i < nReplicas; i++) {
      searchFile = new File(
              parent.getAbsolutePath() + File.separator + i + File.separator + fileNameWithExtension);
      if (!searchFile.exists()) {
        return false;
      }
    }
    return true;
  }

  /**
   * Sample. Restart file write-outs are entirely handled by this method based on how many steps are per cycle.
   * User command input is ignored.
   *
   * @param cycles The number of cycles to run.
   * @param nSteps The number of steps per cycle.
   * @param timeStep The time step.
   * @param printInterval The print interval.
   */
  public void sample(int cycles, long nSteps, double timeStep, double printInterval,
                     double trajInterval, int initTitrDynamics) {
    sample(cycles, nSteps, 0, timeStep, printInterval, trajInterval, initTitrDynamics);
  }

  public void sample(int cycles, long titrSteps, long confSteps, double timeStep,
                     double printInterval, double trajInterval, int initDynamics) {
    done = false;
    terminate = false;
    replica.setRestartFrequency(cycles * (titrSteps + confSteps) * replica.dt + 100); // Full control over restarts handled by this class
    extendedSystem.reGuessLambdas();
    replica.setCoordinates(x);
    int startCycle = 0;
    if (initDynamics > 0 && !restart) {
      extendedSystem.reGuessLambdas();
      extendedSystem.setFixedTitrationState(true);
      extendedSystem.setFixedTautomerState(true);

      logger.info(" ");
      logger.info(" ------------------Start of Equilibration Dynamics------------------\n");
      logger.info(" ");

      if (openMM == null) {
        replica.dynamic(initDynamics, timeStep, printInterval, trajInterval, temp, true, dyn);
      } else {
        x = replica.getCoordinates();
        potential.energy(x);
        openMM.setCoordinates(x);

        openMM.dynamic(initDynamics, timeStep, printInterval, trajInterval, temp, true, dyn);

        x = openMM.getCoordinates();
        replica.setCoordinates(x);
      }

      extendedSystem.setFixedTitrationState(false);
      extendedSystem.setFixedTautomerState(false);
      logger.info(extendedSystem.getLambdaList());
      extendedSystem.writeLambdaHistogram(true);
      extendedSystem.copyESVHistogramTo(parametersHis[rank]); // Copy the ESV hist to be empty

      logger.info(" ");
      logger.info(" ------------------End of Equilibration Dynamics------------------\n");
      logger.info(" ");
    } else if (restart) {
      logger.info(" Omitting initialization steps because this is a restart.");
      startCycle = (int) (restartStep / titrSteps) + 1;
      logger.info(" Restarting pH-REX at cycle " + startCycle + " of " + cycles);
    }

    for (int i = startCycle; i < cycles; i++) {
      // Check for termination request.
      if (terminate) {
        done = true;
        break;
      }

      if (openMM != null) {
        if (confSteps < 3) {
          logger.severe(" Increase number of steps per cycle.");
        }
        dynamicsOpenMM(titrSteps, confSteps, timeStep, printInterval, trajInterval);
      } else {
        dynamics(titrSteps, timeStep, printInterval, trajInterval);
      }
      // Set backups in case job is killed at bad time
      if (i == 0 || backupNeeded) {
        replica.writeRestart();
      }
      copyToBackups();
      replica.writeRestart();

      if (i % 100 == 0) {
        extendedSystem.writeLambdaHistogram(true);
      }

      try {
        logger.info("\n Predicting pKa...");
        logger.info(pkaPred(parametersHis));
      } catch (Exception e) {
        logger.info(e.getMessage());
        logger.warning("Failed to predict pKa. Simulation may be too early on.");
      }

      logger.info(" ");
      logger.info(String.format(" ------------------Exchange Cycle %d------------------\n", i + 1));

      if (i != cycles - 1) {
        exchange();
      }
      logger.info(" ");
      logger.info(" Setting rank " + rank + " esv to pH " + pHScale[rank2Ph[rank]]);
      logger.info(" ");
    }
  }

  private String pkaPred(int[][][] parametersHis) {// parametersHis shape explanation - [rank (w/ some pH)][residue][histogram --> len = 100]

    // Calculate collapsed fractions of 10x10 histograms
    double[][][] collapsedSum = new double[parametersHis.length][parametersHis[0].length][10];
    double[][] collapsedRatio = new double[parametersHis.length][parametersHis[0].length];
    for(int i = 0; i < parametersHis.length; i++) // pH's
    {
      for(int j = 0; j < parametersHis[0].length; j++) // residues
      {
        for(int k = 0; k < 10; k++) // 10 titration states
        {
          double sum = 0;
          for(int l = 0; l < 10; l++) // 10 tautomer states
          {
            sum += parametersHis[i][j][k*10 + l];
          }
          collapsedSum[i][j][k] = sum;
        }
        // Calculate ratio over 0 and 9
        double ratio = collapsedSum[i][j][0] / (collapsedSum[i][j][0] + collapsedSum[i][j][9]);
        collapsedRatio[i][j] = ratio;
      }
    }

    // Calculate hill coeff and pKa's of each residue and print
    StringBuilder output = new StringBuilder();
    double[] n = new double[parametersHis[0].length];
    double[] pka = new double[parametersHis[0].length];
    double[][] residueRatios = new double[parametersHis[0].length][parametersHis.length];
    for(int i = 0; i < parametersHis[0].length; i++) // residues
    {
      // Create array of fractions for this residue
      for(int j = 0; j < parametersHis.length; j++) // pH's
      {
        residueRatios[i][j] = collapsedRatio[j][i];
      }
      //Arrays.sort(residueRatios[i]);

      // L-BFGS minimization of the Henderson-Hasselbalch equation to find the best fit hill coeff and pKa
      double[] temp = TitrationUtils.predictHillCoeffandPka(pHScale, residueRatios[i]);
      n[i] = temp[0];
      pka[i] = temp[1];

      // Print results for this residue
      String residueName = extendedSystem.getTitratingResidueList().get(i).toString();
      output.append(" Residue: ").append(residueName).append("\n");
      output.append(" Fractions (Dep / (Dep + Pro)): ").append(Arrays.toString(residueRatios[i])).append("\n");
      output.append(" pH window: ").append(Arrays.toString(pHScale)).append("\n");
      output.append(" n: ").append(n[i]).append("\n");
      output.append(" pKa: ").append(pka[i]).append("\n");
      output.append("\n");
    }

    return output.toString();
  }

  /**
   * Sets an even pH ladder based on the pH gap.
   */
  public static double[] setEvenSpacePhLadder(double pHGap, double pH, int nReplicas) {
    double range = nReplicas * pHGap;
    double pHMin = pH - range / 2;

    if (nReplicas % 2 != 0) {
      pHMin += pHGap / 2; // Center range if odd num windows
    }

    double[] pHScale = new double[nReplicas];
    for (int i = 0; i < nReplicas; i++) {
      pHScale[i] = pHMin + i * pHGap;
    }

    return pHScale;
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

  /**
   * Evaluate whether to exchange.
   *
   * @param pH The pH to have as the replica target.
   */
  private void compareTwo(int pH) {
    // Ranks for pH A and B
    int rankA;
    int rankB;

    // Load pH, beta and energy for each rank.
    double pHA;
    double pHB;

    double beta;
    double acidostatA;
    double acidostatB;
    double acidostatAatB;
    double acidostatBatA;

    rankA = pH2Rank[pH];
    rankB = pH2Rank[pH + 1];
    int acidAIndex = rank2Ph[rankA] + 1;
    int acidBIndex = rank2Ph[rankB] + 1;

    pHA = parameters[rankA][0];
    pHB = parameters[rankB][0];

    beta = KCAL_TO_GRAM_ANG2_PER_PS2 / (temp * kB);

    acidostatA = parameters[rankA][acidAIndex];
    acidostatB = parameters[rankB][acidBIndex];

    acidostatAatB = parameters[rankA][acidBIndex]; // acidostat of rankA evaluated at the pH of rankB
    acidostatBatA = parameters[rankB][acidAIndex];

    logger.info(" ");
    logger.info(
            " From rank " + rank + ": Comparing ranks " + rankA + " (pH = " + pHA + ") & " + rankB
                    + " (pH = " + pHB + ")");

    // Compute the change in energy over kT (E/kT) for the Metropolis criteria.
    double deltaE = beta * ((acidostatAatB + acidostatBatA) - (acidostatA + acidostatB));
    logger.info("pHA = " + pHA + " " + acidostatA + " " + acidostatAatB);
    logger.info("pHB = " + pHB + " " + acidostatB + " " + acidostatBatA);
    logger.info("exp(" + beta + " * ((" + acidostatAatB + " + " + acidostatBatA + ") - (" + acidostatA + " + " + acidostatB + ")))");
    logger.info(" DeltaE: " + deltaE);

    //Count the number of trials for each temp
    pHTrialCount[pH]++;

    // If the Metropolis criteria is satisfied, do the switch.
    if (deltaE < 0.0 || random.nextDouble() < exp(-deltaE)) {
      int[][] tempHis = new int[extendedSystem.getTitratingResidueList().size()][100];
      for (int i = 0; i < extendedSystem.getTitratingResidueList().size(); i++) {
        System.arraycopy(parametersHis[rankA][i], 0, tempHis[i], 0, tempHis[i].length);
        System.arraycopy(parametersHis[rankB][i], 0, parametersHis[rankA][i], 0,
                parametersHis[rankA][i].length);
        System.arraycopy(tempHis[i], 0, parametersHis[rankB][i], 0, parametersHis[rankB][i].length);
      }
      parameters[rankA][0] = pHB;
      parameters[rankB][0] = pHA;

      //Log the new parameters
      logger.info(" pH After swap " + parameters[rankA][0] + " rankA " + rankA + " acidostat: " + Arrays.toString(parameters[rankA]));
      logger.info(" pH After swap " + parameters[rankB][0] + " rankB " + rankB + " acidostat: " + Arrays.toString(parameters[rankB]));

      // Map pH to process ranks.
      pH2Rank[pH] = rankB;
      pH2Rank[pH + 1] = rankA;

      // Map ranks to pH.
      rank2Ph[rankA] = pH + 1;
      rank2Ph[rankB] = pH;

      pHAcceptedCount[pH]++;
    }
  }

  /** All processes complete the exchanges identically given the same Random number seed. */
  private void exchange() {
    // Loop over top and bottom parts of pH scale
    for (int pH = 0; pH < nReplicas - 1; pH++) {
      compareTwo(pH);
    }

    logger.info(" ");
    // Print Exchange Info
    for (int i = 0; i < pHScale.length - 1; i++) {
      double pHAcceptance = pHAcceptedCount[i] * 100.0 / (pHTrialCount[i]);
      logger.info(
              " Acceptance for pH " + pHScale[i] + " to be exchanged with pH " + pHScale[i + 1] + ": "
                      + pHAcceptance);
    }
    logger.info(" ");
  }

  /**
   * Blocking dynamic steps: when this method returns each replica has completed the requested number
   * of steps. Both OpenMM and CPU implementations exist
   *
   * @param titrSteps the number of time steps on CPU.
   * @param confSteps the number of time steps on GPU.
   * @param timeStep the time step.
   * @param printInterval the number of steps between logging updates.
   * @param saveInterval the number of steps between saving snapshots.
   */
  private void dynamicsOpenMM(final long titrSteps, final long confSteps, final double timeStep,
                              final double printInterval, final double saveInterval) {

    int i = rank2Ph[rank];
    // Update this ranks' parameter array to be consistent with the dynamics.
    myParameters[0] = pHScale[i];
    extendedSystem.setConstantPh(myParameters[0]);
    extendedSystem.copyESVHistogramTo(parametersHis[rank]);

    // Start this processes MolecularDynamics instance sampling.
    boolean initVelocities = true;

    x = replica.getCoordinates();
    potential.energy(x);
    openMM.setCoordinates(x);

    openMM.dynamic(confSteps, timeStep, printInterval, saveInterval, temp, initVelocities, dyn);

    x = openMM.getCoordinates();
    replica.setCoordinates(x);

    double forceWriteInterval = titrSteps * .001;
    replica.dynamic(titrSteps, timeStep, printInterval, forceWriteInterval, temp, initVelocities, dyn);

    reEvaulateAcidostats();

    extendedSystem.getESVHistogram(parametersHis[rank]);

    // Gather all parameters from the other processes.
    try {
      world.allGather(myParametersBuf, parametersBuf);
      world.allGather(parametersHisBuf[rank], parametersHisBuf);
    } catch (IOException ex) {
      String message = " Replica Exchange allGather failed.";
      logger.log(Level.SEVERE, message, ex);
    }
  }


  private void dynamics(long nSteps, double timeStep, double printInterval, double saveInterval) {
    int i = rank2Ph[rank];

    // Update this ranks' parameter array to be consistent with the dynamics.
    myParameters[0] = pHScale[i];
    extendedSystem.setConstantPh(myParameters[0]);

    extendedSystem.copyESVHistogramTo(parametersHis[rank]);

    // Start this processes MolecularDynamics instance sampling.
    boolean initVelocities = true;

    replica.dynamic(nSteps, timeStep, printInterval, saveInterval, temp, initVelocities, dyn);

    reEvaulateAcidostats();

    extendedSystem.getESVHistogram(parametersHis[rank]);

    // Gather all parameters from the other processes.
    try {
      world.allGather(myParametersBuf, parametersBuf);
      world.allGather(parametersHisBuf[rank], parametersHisBuf);
    } catch (IOException ex) {
      String message = " Replica Exchange allGather failed.";
      logger.log(Level.SEVERE, message, ex);
    }
  }

  private void reEvaulateAcidostats() {
    for(int j = 0; j < pHScale.length; j++){
      extendedSystem.setConstantPh(pHScale[j]); // A at A+gap(B)
      myParameters[1+j] = extendedSystem.getBiasEnergy();
    }

    extendedSystem.setConstantPh(myParameters[0]);
  }

  private void copyToBackups() {
    try {
      Files.move(esv.toPath(), esv.toPath().resolveSibling(esvBackup.getName()),
              StandardCopyOption.REPLACE_EXISTING);
    } catch (IOException e) {
      String message = " Could not copy ESV histogram to backup - dynamics terminated.";
      logger.log(Level.WARNING, message);
      throw new RuntimeException(e);
    }

    try {
      Files.move(dyn.toPath(), dyn.toPath().resolveSibling(dynBackup.getName()),
              StandardCopyOption.REPLACE_EXISTING);
    } catch (IOException e) {
      String message = " Could not copy dyn restart to backup - dynamics terminated.";
      logger.log(Level.WARNING, message);
      throw new RuntimeException(e);
    }
  }
}