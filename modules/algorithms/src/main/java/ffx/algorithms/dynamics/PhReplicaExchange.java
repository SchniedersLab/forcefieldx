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
import edu.rit.mp.IntegerBuf;
import edu.rit.mp.buf.IntegerMatrixBuf_1;
import edu.rit.pj.Comm;
import ffx.algorithms.Terminatable;
import ffx.numerics.Potential;
import ffx.potential.extended.ExtendedSystem;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.util.FastMath;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.kB;
import static org.apache.commons.math3.util.FastMath.exp;

/**
 * The ReplicaExchange implements pH replica exchange methods.
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 * @since 1.0
 */
public class PhReplicaExchange implements Terminatable {

  // TODO: Because the program now writes out histograms to restarts after each dynamics, there is no need to communicate the histograms, since they are written out on disk
  // TODO: This means that if we save the map that tells us what rank used to have the pH we are at, we can just read that file back in

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
  private final int[][][] parametersHis;
  /**
   * Each parameter array is wrapped inside a Parallel Java DoubleBuf for the All-Gather
   * communication calls.
   */
  private final DoubleBuf[] parametersBuf;
  private final IntegerBuf[] parametersHisBuf;
  private final MolecularDynamics replica;
  private boolean done = false;
  private boolean restart = true;
  private boolean terminate = false;
  private final double[] myParameters;
  private final DoubleBuf myParametersBuf;
  private final int[] pH2Rank;
  private final int[] rank2Ph;
  private final double[] pHScale;
  private final int[] pHAcceptedCount;
  private final int[] rankAcceptedCount;
  private final int[] pHTrialCount;
  private final double temp;
  private final ExtendedSystem extendedSystem;
  private final double pH;
  private final double gapSize;
  private double[] x;
  private MolecularDynamicsOpenMM openMM = null;
  private Potential potential;
  private File dyn;
  private int restartStep;

  /**
   * ReplicaExchange constructor.
   *
   * @param molecularDynamics a {@link MolecularDynamics} object.
   * @param pH pH = pKa <-- will be changed from this initial value
   * @param extendedSystem extended system attached to this process
   * @param pHGap the gap in pH units between replicas
   * @param temp temperature of replica
   */
  public PhReplicaExchange(
      MolecularDynamics molecularDynamics, File structureFile, double pH, double pHGap, double temp, ExtendedSystem extendedSystem) {

    this.replica = molecularDynamics;
    this.temp = temp;
    this.extendedSystem = extendedSystem;
    this.pH = pH;
    this.gapSize = pHGap;

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

    setEvenSpacePhLadder(pHGap);
    extendedSystem.setConstantPh(pHScale[rank]);

    // File stuff
    File rankDir = new File(structureFile.getParent() + File.separator + rank);
    restart = !rankDir.mkdir();

    File esv = new File(rankDir.getPath() + File.separator + FilenameUtils.removeExtension(structureFile.getName()) + ".esv");
    dyn = new File(rankDir.getPath() + File.separator + FilenameUtils.removeExtension((structureFile.getName())) + ".dyn");
    extendedSystem.setESVFile(esv);
    molecularDynamics.setFallbackDynFile(dyn);

    restartStep = 0;
    // Count how far the restart should start at and make sure that each dir has the correct files and set up maps
    if(restart){
      ArrayList<Double> readPhScale = new ArrayList<>();
      for(int i = 0; i < nReplicas; i++){
        File checkESV = new File(rankDir.getParent() + File.separator + i + File.separator + FilenameUtils.removeExtension(structureFile.getName()) + ".esv");
        File checkDYN = new File(rankDir.getParent() + File.separator + i + File.separator + FilenameUtils.removeExtension(structureFile.getName()) + ".dyn");
        if(!(checkDYN.exists() && checkESV.exists())) {
          logger.severe("Rank " + i + " does not contain the correct .dyn or .esv.");
        }

        try(BufferedReader br = new BufferedReader(new FileReader(checkESV))){
          String data = br.readLine();
          //TODO: Think about how this works with protein restarts
          while(data != null){
            List<String> tokens = Arrays.asList(data.split(" +"));
            if(tokens.contains("pH:")){
              double pHOfRankI = Double.parseDouble(tokens.get(tokens.indexOf("pH:") + 1));
              if(readPhScale.contains(pHOfRankI)){
                logger.warning(" Duplicate pH values found. ");
                logger.warning(" Restart is unable to run due to a loss of data on program killing. ");
                restart = false;
                break;
              }
              else {
                readPhScale.add(pHOfRankI);
              }
              // Set up map
              for(int j = 0; j < pHScale.length; j++){
                if(pHScale[j] == pHOfRankI){
                  rank2Ph[i] = j;
                  pH2Rank[j] = i;
                }
              }

              br.readLine();
              br.readLine();
              int sum = 0;
              for(int j = 0; j < 10; j++){ // 10 titr windows
                data = br.readLine().trim();
                tokens = List.of(data.split(" +"));
                for(int k = 0; k < 10; k++){ // 10 tautomer windows
                  sum += Integer.parseInt(tokens.get(k + 1));
                }
              }
              if(i == 0) {
                restartStep = sum;
                logger.info(" Restart already completed " + restartStep + " steps."); // Restart step loses one every restart???
              } else if(restartStep != sum){
                logger.warning(" Restart received uneven sums. Starting from the lowest one. Some windows may have more data then expected");
                logger.info("Restart Step Current: " + restartStep);
                restartStep = Math.min(sum, restartStep);
                logger.info(" Restart Step New: " + restartStep);
              }
            }
            data = br.readLine();
          }
        } catch (IOException e) {
          e.printStackTrace();
          logger.severe(" Error reading " + checkESV.getAbsolutePath());
        }
        catch (IndexOutOfBoundsException e){
          e.printStackTrace();
          logger.warning(" pH values changed or could not be found. Running from start");
          restart = false;
        }
        if(!restart){
          logger.severe(" Was not able to restart...");
        }
      }
    }


    // Set restart files
    if(restart) {
      logger.info(" Rank " + rank + " staring at pH " + pHScale[rank2Ph[rank]]);
      extendedSystem.setConstantPh(pHScale[rank2Ph[rank]]);
      extendedSystem.setESVFile(esv);
      molecularDynamics.setFallbackDynFile(dyn);
    }

    random = new Random();
    random.setSeed(0);

    // Create arrays to store the parameters of all processes.
    parameters = new double[nReplicas][4];
    parametersHis = new int[nReplicas][extendedSystem.getTitratingResidueList().size()][100];
    parametersBuf = new DoubleBuf[nReplicas];
    parametersHisBuf = new IntegerBuf[nReplicas];
    for (int i = 0; i < nReplicas; i++) {
      parametersBuf[i] = DoubleBuf.buffer(parameters[i]);
      parametersHisBuf[i] = IntegerMatrixBuf_1.buffer(parametersHis[i]);
    }

    // Set the ESV that this rank will share has correct values after a restart
    if(restart){
      extendedSystem.getESVHistogram(parametersHis[rank]);
    }

    // A convenience reference to the parameters of this process are updated
    // during communication calls.
    myParameters = parameters[rank];
    //myParametersHis = parametersHis[rank];
    myParametersBuf = parametersBuf[rank];
    //myParametersHisBuf = parametersHisBuf[rank];
  }

  /**
   * OpenMM ReplicaExchange constructor.
   *
   * @param molecularDynamics a {@link MolecularDynamics} object.
   * @param pH pH = pKa <-- will be changed from this initial value
   * @param extendedSystem extended system attached to this process
   * @param pHGap the gap in pH units between replicas
   * @param temp temperature of replica
   * @param x array of coordinates
   * @param molecularDynamicsOpenMM for running config steps on GPU
   */
  public PhReplicaExchange(
          MolecularDynamics molecularDynamics, File structureFile, double pH, double pHGap, double temp, ExtendedSystem extendedSystem,
          double[] x, MolecularDynamicsOpenMM molecularDynamicsOpenMM, Potential potential) {

    this.replica = molecularDynamics;
    this.temp = temp;
    this.extendedSystem = extendedSystem;
    this.pH = pH;
    this.gapSize = pHGap;
    this.x = x;
    this.openMM = molecularDynamicsOpenMM;
    this.potential = potential;

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

    setEvenSpacePhLadder(pHGap);
    extendedSystem.setConstantPh(pHScale[rank]);

    // File stuff
    File rankDir = new File(structureFile.getParent() + File.separator + rank);
    restart = !rankDir.mkdir();

    File esv = new File(rankDir.getPath() + File.separator + FilenameUtils.removeExtension(structureFile.getName()) + ".esv");
    dyn = new File(rankDir.getPath() + File.separator + FilenameUtils.removeExtension((structureFile.getName())) + ".dyn");
    extendedSystem.setESVFile(esv);
    molecularDynamics.setFallbackDynFile(dyn);
    molecularDynamicsOpenMM.setAutomaticWriteouts(false);

    restartStep = 0;
    // Count how far the restart should start at and make sure that each dir has the correct files and set up maps
    if(restart){
      ArrayList<Double> readPhScale = new ArrayList<>();
      for(int i = 0; i < nReplicas; i++){
        File checkESV = new File(rankDir.getParent() + File.separator + i + File.separator + FilenameUtils.removeExtension(structureFile.getName()) + ".esv");
        File checkDYN = new File(rankDir.getParent() + File.separator + i + File.separator + FilenameUtils.removeExtension(structureFile.getName()) + ".dyn");
        if(!(checkDYN.exists() && checkESV.exists())) {
          logger.severe("Rank " + i + " does not contain the correct .dyn or .esv.");
        }

        try(BufferedReader br = new BufferedReader(new FileReader(checkESV))){
          String data = br.readLine();
          while(data != null){
            List<String> tokens = Arrays.asList(data.split(" +"));
            if(tokens.contains("pH:")){
              double pHOfRankI = Double.parseDouble(tokens.get(tokens.indexOf("pH:") + 1));
              if(readPhScale.contains(pHOfRankI)){
                logger.warning(" Duplicate pH values found. ");
                logger.warning(" Restart is unable to run due to a loss of data on program killing. ");
                restart = false;
                break;
              }
              else {
                readPhScale.add(pHOfRankI);
              }
              // Set up map
              for(int j = 0; j < pHScale.length; j++){
                if(pHScale[j] == pHOfRankI){
                  rank2Ph[i] = j;
                  pH2Rank[j] = i;
                }
              }

              br.readLine();
              br.readLine();
              int sum = 0;
              for(int j = 0; j < 10; j++){ // 10 titr windows
                data = br.readLine().trim();
                tokens = List.of(data.split(" +"));
                for(int k = 0; k < 10; k++){ // 10 tautomer windows
                  sum += Integer.parseInt(tokens.get(k + 1));
                }
              }
              if(i == 0) {
                restartStep = sum;
                logger.info(" Restart already completed " + restartStep + " steps."); // Restart step loses one every restart???
              } else if(restartStep != sum){
                logger.warning(" Restart received uneven sums. Starting from the lowest one. Some windows may have more data then expected");
                logger.info("Restart Step Current: " + restartStep);
                restartStep = Math.min(sum, restartStep);
                logger.info(" Restart Step New: " + restartStep);
              }
            }
            data = br.readLine();
          }
        } catch (IOException e) {
          e.printStackTrace();
          logger.severe(" Error reading " + checkESV.getAbsolutePath());
        }
        catch (IndexOutOfBoundsException e){
          e.printStackTrace();
          logger.warning(" pH values changed or could not be found. Running from start");
          restart = false;
        }
        if(!restart){
          logger.severe(" Was not able to restart...");
        }
      }
    }


    // Set restart files
    if(restart) {
      logger.info(" Rank " + rank + " staring at pH " + pHScale[rank2Ph[rank]]);
      extendedSystem.setConstantPh(pHScale[rank2Ph[rank]]);
      extendedSystem.setESVFile(esv);
      logger.info("Rank " + rank + " starting hist:");
      extendedSystem.writeLambdaHistogram();
      molecularDynamics.setFallbackDynFile(dyn);
    }

    random = new Random();
    random.setSeed(0);

    // Create arrays to store the parameters of all processes.
    parameters = new double[nReplicas][4];
    parametersHis = new int[nReplicas][extendedSystem.getTitratingResidueList().size()][100];
    parametersBuf = new DoubleBuf[nReplicas];
    parametersHisBuf = new IntegerMatrixBuf_1[nReplicas];
    for (int i = 0; i < nReplicas; i++) {
      parametersBuf[i] = DoubleBuf.buffer(parameters[i]);
      parametersHisBuf[i] = IntegerMatrixBuf_1.buffer(parametersHis[i]);
    }

    // Set the ESV that this rank will share has correct values after a restart
    if(restart){
      extendedSystem.getESVHistogram(parametersHis[rank]);
    }

    // A convenience reference to the parameters of this process are updated
    // during communication calls.
    myParameters = parameters[rank];
    //myParametersHis = parametersHis[rank];
    myParametersBuf = parametersBuf[rank];
    //myParametersHisBuf = parametersHisBuf[rank];
  }

  /**
   * Sample.
   *
   * @param cycles        a int.
   * @param nSteps        a int.
   * @param timeStep      a double.
   * @param printInterval a double.
   */
  public void sample(int cycles, long nSteps, double timeStep, double printInterval, double trajInterval, int initTitrDynamics) {
    sample(cycles, nSteps, 0, timeStep, printInterval, trajInterval, initTitrDynamics);
  }

  public void sample(int cycles, long titrSteps, long confSteps, double timeStep, double printInterval, double trajInterval, int initTitrDynamics) {
    done = false;
    terminate = false;
    replica.setRestartFrequency(cycles * (titrSteps + confSteps) * replica.dt + 100);

    if(extendedSystem.guessTitrState){
      extendedSystem.reGuessLambdas();
    }

    if(initTitrDynamics > 0 && !restart) {
      replica.dynamic(initTitrDynamics, timeStep,
              printInterval, trajInterval, temp, true, dyn);
      logger.info(extendedSystem.getLambdaList());
      extendedSystem.writeLambdaHistogram();
      extendedSystem.copyESVHistogramTo(parametersHis[rank]); // Copy the ESV hist to be empty

      logger.info(" ");
      logger.info(" ------------------End of Equilibration Dynamics------------------\n");
      logger.info(" ");
    }

    int startCycle = 0;
    if(restart){
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

      if(openMM != null){
        if(confSteps < 3)
        {
          logger.severe("Increase number of steps per cycle.");
        }
        dynamicsOpenMM(titrSteps + 2, confSteps, timeStep, printInterval, trajInterval);
      }
      else {
        dynamics(titrSteps + 1, timeStep, printInterval, trajInterval);
      }
      replica.writeRestart();

      logger.info(" ");
      logger.info(String.format(" ------------------Exchange Cycle %d------------------\n", i+1));

      if(i != cycles-1) {
        exchange();
      }
      logger.info(" ");
      logger.info(" Setting rank " + rank + " esv to pH " + pHScale[rank2Ph[rank]]);
      logger.info(" ");
    }
  }

  public double[] getpHScale(){
    return pHScale;
  }

  /**
   * setEvenPhLadder.
   *
   */
  private void setEvenSpacePhLadder(double pHGap){
    double range = world.size() * pHGap;
    double pHMin = pH - range/2;

    if(nReplicas % 2 != 0){
      pHMin += pHGap/2;
    }

    for(int i = 0; i < nReplicas; i++){
      pHScale[i] = pHMin + i * pHGap;
      rank2Ph[i] = i;
      pH2Rank[i] = i;
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

  /**
   * Evaluate whether or not to exchange
   * @param pH what pH to have as the replica target
   */
  private void compareTwo(int pH){
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

    pHA = parameters[rankA][0];
    pHB = parameters[rankB][0];

    beta = KCAL_TO_GRAM_ANG2_PER_PS2 / (temp * kB);

    acidostatA = parameters[rankA][2];
    acidostatB = parameters[rankB][2];

    acidostatAatB = parameters[rankA][3]; // acidostat of rankA evaluated at the pH of rankB
    acidostatBatA = parameters[rankB][1];

    logger.info(" ");
    logger.info(" From rank " + rank + ": Comparing ranks " + rankA + " (pH = " + pHA + ") & " + rankB + " (pH = " + pHB + ")");

    // Compute the change in energy over kT (E/kT) for the Metropolis criteria.
    double deltaE = beta * ((acidostatAatB + acidostatBatA) - (acidostatA + acidostatB));
    logger.info(" deltaE = " + beta + " * ((" + acidostatAatB + " + " + acidostatBatA + ") - (" + acidostatA + " + " + acidostatB + "))");
    logger.info(" DeltaE: " + deltaE);

    //Count the number of trials for each temp
    pHTrialCount[pH]++;

    // If the Metropolis criteria is satisfied, do the switch.
    if (deltaE < 0.0 || random.nextDouble() < exp(-deltaE)) {
      int[][] tempHis = new int[extendedSystem.getTitratingResidueList().size()][100];
      for (int i = 0; i < extendedSystem.getTitratingResidueList().size(); i++) {
        System.arraycopy(parametersHis[rankA][i],0,tempHis[i],0,tempHis[i].length);
        System.arraycopy(parametersHis[rankB][i],0,parametersHis[rankA][i],0,parametersHis[rankA][i].length);
        System.arraycopy(tempHis[i],0,parametersHis[rankB][i],0,parametersHis[rankB][i].length);
      }
      // Swap pH and energy values.
      parameters[rankA][0] = pHB;
      parameters[rankB][0] = pHA;

      // Map pH to process ranks.
      pH2Rank[pH] = rankB;
      pH2Rank[pH + 1] = rankA;

      // Map ranks to pH.
      rank2Ph[rankA] = pH + 1;
      rank2Ph[rankB] = pH;

      pHAcceptedCount[pH]++;
      if (pHA < pHB){
        rankAcceptedCount[rankA]++;
      } else {
        rankAcceptedCount[rankB]++;
      }
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
      logger.info(" Acceptance for pH " + pHScale[i] + " to be exchanged with pH " + pHScale[i+1] + ": " + pHAcceptance);
      }
    logger.info(" ");
  }

  /**
   * Blocking dynamic steps: when this method returns each replica has completed the requested
   * number of steps. Both OpenMM and CPU implementations exist
   *
   * @param titrSteps the number of time steps on CPU.
   * @param confSteps the number of time steps on GPU.
   * @param timeStep the time step.
   * @param printInterval the number of steps between loggging updates.
   * @param saveInterval the number of steps between saving snapshots.
   */
  private void dynamicsOpenMM(final long titrSteps, final long confSteps, final double timeStep, final double printInterval, final double saveInterval) {

    int i = rank2Ph[rank];
    extendedSystem.setConstantPh(pHScale[i]);
    extendedSystem.copyESVHistogramTo(parametersHis[rank]);

    // Start this processes MolecularDynamics instance sampling.
    boolean initVelocities = true;

    int titrStepsOne = (int) titrSteps / 2;
    int titrStepsTwo = (int) FastMath.ceil(titrSteps / 2.0);

    replica.dynamic(titrStepsOne, timeStep, printInterval, saveInterval, temp, initVelocities, dyn);

    x = replica.getCoordinates();
    potential.energy(x);
    openMM.setCoordinates(x);

    openMM.dynamic(confSteps, timeStep, printInterval, saveInterval, temp, initVelocities, dyn);

    x = openMM.getCoordinates();
    replica.setCoordinates(x);

    replica.dynamic(titrStepsTwo, timeStep, printInterval, saveInterval, temp, initVelocities, dyn);

    // Update this ranks' parameter array to be consistent with the dynamics.

    myParameters[0] = pHScale[i];
    myParameters[2] = extendedSystem.getBiasEnergy();
    logger.info(" ");

    // Evaluate acidostat of ES at different pHs
    extendedSystem.setConstantPh(myParameters[0] - gapSize);
    myParameters[1] = extendedSystem.getBiasEnergy();

    extendedSystem.setConstantPh(myParameters[0] + gapSize);
    myParameters[3] = extendedSystem.getBiasEnergy();

    extendedSystem.setConstantPh(myParameters[0]);

    extendedSystem.getESVHistogram(parametersHis[rank]);
    extendedSystem.writeLambdaHistogram();

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

    extendedSystem.setConstantPh(pHScale[i]);

    extendedSystem.copyESVHistogramTo(parametersHis[rank]);

    // Start this processes MolecularDynamics instance sampling.
    boolean initVelocities = true;

    replica.dynamic(
            nSteps, timeStep, printInterval, saveInterval, temp, initVelocities, dyn);

    // Update this ranks' parameter array to be consistent with the dynamics.
    myParameters[0] = pHScale[i];
    myParameters[2] = extendedSystem.getBiasEnergy();

    // Evaluate acidostat of ES at different pHs
    extendedSystem.setConstantPh(myParameters[0] - gapSize); // B at A-gap
    myParameters[1] = extendedSystem.getBiasEnergy();

    extendedSystem.setConstantPh(myParameters[0] + gapSize); // A at A+gap(B)
    myParameters[3] = extendedSystem.getBiasEnergy();

    extendedSystem.setConstantPh(myParameters[0]);

    extendedSystem.getESVHistogram(parametersHis[rank]);
    extendedSystem.writeLambdaHistogram();

    // Gather all parameters from the other processes.
    try {
      world.allGather(myParametersBuf, parametersBuf);
      world.allGather(parametersHisBuf[rank], parametersHisBuf);
    } catch (IOException ex) {
      String message = " Replica Exchange allGather failed.";
      logger.log(Level.SEVERE, message, ex);
    }
  }
}