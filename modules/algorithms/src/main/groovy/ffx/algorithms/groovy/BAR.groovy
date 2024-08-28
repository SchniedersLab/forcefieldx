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
package ffx.algorithms.groovy

import edu.rit.mp.DoubleBuf
import edu.rit.mp.IntegerBuf
import edu.rit.pj.Comm
import edu.rit.pj.ParallelTeam
import ffx.algorithms.cli.AlgorithmsScript
import ffx.crystal.Crystal
import ffx.crystal.CrystalPotential
import ffx.numerics.Potential
import ffx.numerics.estimator.BennettAcceptanceRatio
import ffx.numerics.estimator.EstimateBootstrapper
import ffx.numerics.estimator.SequentialEstimator
import ffx.numerics.math.BootStrapStatistics
import ffx.potential.MolecularAssembly
import ffx.potential.Utilities
import ffx.potential.bonded.LambdaInterface
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.TopologyOptions
import ffx.potential.parsers.BARFilter
import ffx.potential.parsers.SystemFilter
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.configuration2.Configuration
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.logging.Level

import static ffx.utilities.Constants.NS2SEC
import static java.lang.Integer.MIN_VALUE
import static java.lang.Integer.parseInt;
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.min

/**
 * The BAR script find the free energy difference across a lambda window. It presently assumes
 * that the number of files composing the first end of the window equals the number of files
 * composing the other end.
 * <br>
 * Usage:
 * <br>
 * ffxc BAR [options] &lt;structures1&gt &lt;structures2&gt;
 */
@Command(description = " Evaluates a free energy change with the Bennett Acceptance Ratio algorithm using pregenerated snapshots.", name = "BAR")
class BAR extends AlgorithmsScript {

  @Mixin
  private AlchemicalOptions alchemical

  @Mixin
  private TopologyOptions topology

  @Option(names = ["--l2", "--lambdaTwo"], paramLabel = "1.0",
          description = "Lambda value for the upper edge of the window")
  private double lambda2 = 1.0

  @Option(names = ["-t", "--temperature"], paramLabel = "298.15",
          description = "Temperature for system")
  private double temperature = 298.15

  @Option(names = ["--dV", "--volume"], paramLabel = "false",
          description = "Write out snapshot volumes to the Tinker BAR file.")
  private boolean includeVolume = false

  @Option(names = ["--tb", "--tinkerBAR"], paramLabel = "false",
          description = "Write out a Tinker BAR file.")
  private boolean tinkerBAR = false

  @Option(names = ["--nw", "--nWindows"], paramLabel = "-1",
          description = "If set, auto-determine lambda values and subdirectories (overrides other flags).")
  private int nWindows = -1

  @Option(names = ["--utb", "--useTinker"], paramLabel = "false",
          description = "If set, use tinker BAR files for energy snapshots.")
  private boolean useTinkerBAR = false

  @Option(names = ["--sa", "--sortedArc"], paramLabel = "false",
          description = "If set, use sorted archive values.")
  private boolean sortedArc = false

  @Option(names = ["--ss", "--startSnapshot"], paramLabel = "0",
          description = "Start at this snapshot when reading in tinker BAR files.")
  private int startingSnapshot = 0

  @Option(names = ["--es", "--endSnapshot"], paramLabel = "0",
          description = "End at this snapshot when reading in tinker BAR files.")
  private int endingSnapshot = 0

  /**
   * --ni or --nIterattions Maximum number of allowable iterations for BAR calculation.
   */
  @Option(names = ["--ni", "--nIterations"], paramLabel = "100",
          description = "Specify the maximum number of iterations for BAR convergence.")
  private int nIterations = 100

  /**
   * -e or --eps Convergence criterion for BAR iteration..
   */
  @Option(names = ["-e", "--eps"], paramLabel = "1.0E-7",
          description = "Specify convergence cutoff for BAR calculation.")
  private double eps = 1.0E-7

  @Option(names = ["--lambdaArray"], paramLabel = "0,0.2,0.25,...,1.0",
          description = "Custom lambda values as a comma-separated list.")
  private String lambdaArrayStr;

  private double[] lambdaArray;

  public void parseLambdaArray() {
    if (lambdaArrayStr != null && !lambdaArrayStr.isEmpty()) {
      String[] tokens = lambdaArrayStr.split(",");
      lambdaArray = new double[tokens.length];
      for (int i = 0; i < tokens.length; i++) {
        lambdaArray[i] = Double.parseDouble(tokens[i].trim());
      }
    }
  }

  /**
   * The final argument(s) should be filenames for lambda windows in order.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
          description = 'A single PDB/XYZ when windows are auto-determined (or two for dual topology). Two trajectory files for BAR between two ensembles (or four for dual topology).')
  List<String> filenames = null

  /**
   * Number of threads available to utilize.
   */
  private int threadsAvail = ParallelTeam.getDefaultThreadCount()
  /**
   * Threads per walker
   */
  private int threadsPer = threadsAvail
  /**
   * Molecular assemblies to compute energies.
   */
  MolecularAssembly[] molecularAssemblies
  /**
   * Potential object for the crystal
   */
  CrystalPotential potential
  /**
   * File openers to read system information from files.
   */
  SystemFilter[] openers
  /**
   * Openers for existing BAR files.
   */
  BARFilter[] barOpeners
  /**
   * Writers to create BAR files.
   */
  BARFilter[] barWriters
  /**
   * Location of windows
   */
  List<String> windowDirectories = new ArrayList<>()
  /**
   * Additional properties
   */
  private Configuration additionalProperties
  /**
   * Input files for BAR.
   */
  private String[] files
  /**
   * Free energy difference from forward FEP.
   */
  private double forwardFEP
  /**
   * Enthalpy difference from forward FEP.
   */
  private double forwardEnthalpy
  /**
   * Entropy difference from forward FEP.
   */
  private double forwardEntropy
  /**
   * Free energy difference from backward FEP.
   */
  private double backwardFEP
  /**
   * Enthalpy difference from backward FEP.
   */
  private double backwardEnthalpy
  /**
   * Entropy difference from backward FEP.
   */
  private double backwardEntropy
  /**
   * Free energy difference from BAR.
   */
  private double barEnergy
  /**
   * Free energy difference from BAR using boot-strapping.
   */
  private double barEnergyBS
  /**
   * Enthalpy difference from BAR.
   */
  private double barEnthalpy
  /**
   * Entropy difference from BAR.
   */
  private double barEntropy
  /**
   * Number of Topologies
   */
  private int numTopologies
  boolean autodetect = false
  double[] lambdaValues
  int nFiles
  /**
   * Parallel Java world communicator.
   */
  private Comm world;
  /**
   * If false, do not use MPI communication.
   */
  private boolean useMPI;
  /**
   * Number of processes.
   */
  private int numProc;
  /**
   * Rank of this process.
   */
  private int rank;
  /**
   * The energy matrix stores a single energy value from each process. The array is of size
   * [numProc][numWorkItems][snapshots].
   */
  private double[][][] energiesLowPJ;
  /**
   * The energy matrix stores a single energy value from each process. The array is of size
   * [numProc][numWorkItems][snapshots].
   */
  private double[][][] energiesAtPJ;
  /**
   * The energy matrix stores a single energy value from each process. The array is of size
   * [numProc][numWorkItems][snapshots].
   */
  private double[][][] energiesHighPJ;
  /**
   * The volume matrix stores a single volume value from each process. The array is of size
   * [numProc][numWorkItems][snapshots].
   */
  private double[][][] volumePJ;
  /**
   * The number of models matrix stores a single volume value from each process. The array is of size
   * [numProc][1].
   */
  private int[][] maxModelsPJ
  /**
   * The energy matrix stores a single energy value from each process. The array is of size
   * [nWindows][nSnapshots].
   */
  private double[][] energiesLow;
  /**
   * The energy matrix stores a single energy value from each process. The array is of size
   * [nWindows][nSnapshots].
   */
  private double[][] energiesAt;
  /**
   * The energy matrix stores a single energy value from each process. The array is of size
   * [nWindows][nSnapshots].
   */
  private double[][] energiesHigh;
  /**
   * The volume matrix stores a single volume value from each process. The array is of size
   * [nWindows][nSnapshots].
   */
  private double[][] volume;
  /**
   * Each distance is wrapped inside a DoubleBuf for MPI communication.
   */
  private DoubleBuf[] buffersLow;
  /**
   * Each distance is wrapped inside a DoubleBuf for MPI communication.
   */
  private DoubleBuf[] buffersAt;
  /**
   * Each distance is wrapped inside a DoubleBuf for MPI communication.
   */
  private DoubleBuf[] buffersHigh;
  /**
   * Each distance is wrapped inside a DoubleBuf for MPI communication.
   */
  private DoubleBuf[] buffersVolume;
  /**
   * Each distance is wrapped inside a DoubleBuf for MPI communication.
   */
  private IntegerBuf[] buffersMax;
  /**
   * Convenience reference for the DoubleBuf of this process.
   */
  private DoubleBuf myBufferLow;
  /**
   * Convenience reference for the DoubleBuf of this process.
   */
  private DoubleBuf myBufferAt;
  /**
   * Convenience reference for the DoubleBuf of this process.
   */
  private DoubleBuf myBufferHigh;
  /**
   * Convenience reference for the DoubleBuf of this process.
   */
  private DoubleBuf myBufferVolume;
  /**
   * Convenience reference for the DoubleBuf of this process.
   */
  private IntegerBuf myBufferMax;
  /**
   * The amount of work based on windows for each process.
   */
  private int numWorkItems;

  /**
   * Maximum number of trials to be used for bootstrap.
   */
  final long MAX_BOOTSTRAP_TRIALS = 100000L

  /**
   * Sets an optional Configuration with additional properties.
   * @param additionalProps
   */
  void setProperties(Configuration additionalProps) {
    this.additionalProperties = additionalProps
  }

  /**
   * BAR Constructor.
   */
  BAR() {
    this(new Binding())
  }

  /**
   * BAR Constructor.
   * @param binding The Groovy Binding to use.
   */
  BAR(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  BAR run() {
    // Begin boilerplate code.
    if (!init()) {
      return this
    }

    /**
     * Load user supplied files into an array.
     */
    nFiles = filenames.size()
    files = new String[nFiles]
    for (int i = 0; i < nFiles; i++) {
      files[i] = filenames.get(i)
    }

    numTopologies = 1
    if (lambdaArrayStr != null) {
      parseLambdaArray();
      nWindows=lambdaArray.size();
      System.out.println("Lambda Array: " + Arrays.toString(lambdaArray));
    } else {
      System.out.println("Lambda Array is null");
    }

    if (nFiles <= 1 && nWindows < 2 && lambdaArray == null) {
      logger.info(' At least two ensembles must be specified')
      return this
    } else if (nFiles == 1 && nWindows >= 2) {
      logger.
              info(format(' Auto-detecting %d windows for single topology:\n %s.', nWindows, files[0]))
    } else if (nFiles == 2 && nWindows >= 2) {
      logger.info(format(' Auto-detecting %d windows for dual topology:\n %s\n %s.', nWindows, files[0], files[1]))
      numTopologies = 2
    } else if (nFiles == 2 && nWindows < 2) {
      logger.info(format(' Applying BAR between two single topology ensembles:\n %s\n %s.', files[0], files[1]))
    } else if (nFiles == 4 && nWindows < 2) {
      logger.info(format(' Applying BAR between two dual topology ensembles:\n %s %s\n %s %s.', files[0], files[1], files[2], files[3]))
      numTopologies = 2
    } else {
      logger.info(format(" Inconsistent input of files (%3d) and/or windows (%3d).", nFiles, nWindows))
      return this
    }


    if (lambdaArray != null && lambdaArray.size() > 0) {
      autodetect = true
      nWindows = lambdaArray.size();
      lambdaValues = lambdaArray as double[];
      for (int i = 0; i < nWindows; i++) {
        for (int j = 0; j < nFiles; j++) {
          String fullPathToFile = FilenameUtils.getFullPath(files[j]);
          String directoryFullPath = fullPathToFile.replace(files[j], "") + i;
          windowDirectories.add(directoryFullPath + File.separator + i);
        }
      }
    } else if (nWindows > 1) {
      autodetect = true;
      // Auto-determine subdirectories and their lambda values.
      for (int i = 0; i < nWindows; i++) {
        for (int j = 0; j < nFiles; j++) {
          String fullPathToFile = FilenameUtils.getFullPath(files[j]);
          String directoryFullPath = fullPathToFile.replace(files[j], "") + i;
          windowDirectories.add(directoryFullPath + File.separator + i);
        }
      }
      lambdaValues = new double[nWindows];
      for (int i = 0; i < nWindows; i++) {
        lambdaValues[i] = alchemical.getInitialLambda(nWindows, i, true);
      }
    } else {
      // Default case for when nWindows is not set and no lambdaArray is provided.
      lambdaValues = new double[2];
      lambdaValues[0] = alchemical.getInitialLambda();
      lambdaValues[1] = lambda2;
      nWindows = 2;
    }

    // Could set "getInitialLambda"'s quiet flag to false, but better logging here?
    logger.info(" Lambda values for each window: ")
    int nLambda = lambdaValues.length;
    for (int i = 0; i < nLambda; i++) {
      double l = lambdaValues[i]
      logger.info(format(" Window %3d: %6.4f", i, l))
    }

    // If reading Tinker BAR file(s), allocate storage for openers.
    if (useTinkerBAR) {
      barOpeners = new BARFilter[nWindows]
    }

    // If saving Tinker BAR file(s), allocate storage for writer(s).
    if (tinkerBAR) {
      barWriters = new BARFilter[nWindows]
    }

    // Allocate space for each topology
    molecularAssemblies = new MolecularAssembly[numTopologies]
    openers = new SystemFilter[numTopologies]

    int numParallel = topology.getNumParallel(threadsAvail, numTopologies)
    threadsPer = (int) (threadsAvail / numParallel)

    /*
    Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
    Checking numTopologies == 1 should only be done for scripts that imply
    some sort of lambda scaling.
    The Minimize script, for example, may be running on a single, unscaled physical topology.
    */
    boolean lambdaTerm = (numTopologies == 1 || alchemical.hasSoftcore() || topology.hasSoftcore())

    if (lambdaTerm) {
      System.setProperty("lambdaterm", "true")
    }

    // Relative free energies via the DualTopologyEnergy class require different
    // default free energy parameters than absolute free energies.
    if (numTopologies == 2) {
      // Ligand vapor electrostatics are not calculated. This cancels when the
      // difference between protein and water environments is considered.
      System.setProperty("ligand-vapor-elec", "false")
    }

    if (numTopologies == 2) {
      logger.info(format(" Initializing two topologies for each window."))
    } else {
      logger.info(format(" Initializing a single topology for each window."))
    }

    // Open assemblies
    for (int i = 0; i < numTopologies; i++) {
      MolecularAssembly ma = alchemical.openFile(algorithmFunctions, topology, threadsPer, filenames[i], i)
      molecularAssemblies[i] = ma
      openers[i] = algorithmFunctions.getFilter()
    }

    StringBuilder sb = new StringBuilder(format(
            "\n Using BAR to analyze a free energy change for %s\n ", filenames))
    potential = (CrystalPotential) topology.assemblePotential(molecularAssemblies, threadsAvail, sb)
    Crystal unitCell = potential.getCrystal().getUnitCell()

    boolean isPBC = includeVolume && !unitCell.aperiodic()
    isPBC = isPBC && !unitCell.aperiodic()

    int nSymm = 0
    if (isPBC) {
      nSymm = unitCell.getNumSymOps()
    }


    File file = new File(files[0])
    String directoryPath = file.getAbsoluteFile().getParent() + File.separator
    String[][] fullFilePaths = new String[nWindows][nFiles]
    // Loop over user supplied files.
    for (int j = 0; j < nFiles; j++) {
      // Loop over ensembles.
      for (int i = 0; i < nWindows; i++) {
        String archiveName
        if (sortedArc) {
          archiveName = FilenameUtils.getBaseName(files[j]) + "_E" + i.toString() + ".arc"
        } else {
          archiveName = FilenameUtils.getBaseName(files[j]) + ".arc"
        }
        if (useTinkerBAR && i != nWindows - 1) {
          // Path to Tinker BAR files.
          fullFilePaths[i][j] = directoryPath + "barFiles" + File.separator + "energy_" + i + ".bar"
        } else if (!autodetect) {
          // Path to a file in the same directory as supplied archives.
          fullFilePaths[i][j] = directoryPath + File.separator + archiveName
        } else {
          // Paths to auto-detected subdirectories.
          fullFilePaths[i][j] = directoryPath + i + File.separator + archiveName
        }
      }
    }

    // Load properties file to determine if parallel environment if specified.
    CompositeConfiguration properties = algorithmFunctions.getActiveAssembly().getProperties();
    useMPI = properties.getBoolean("pj.use.mpi", true);
    // Set up parallel objects if necessary.
    if (useMPI) {
      world = Comm.world();
      // Number of processes is equal to world size (often called size).
      numProc = world.size();
      // Each processor gets its own rank (ID of sorts).
      rank = world.rank();

      // Padding of the target array size (inner loop limit) is for parallelization.
      // Target conformations are parallelized over available nodes.
      // For example, if numProc = 8 and nWindows = 12, then paddednWindows = 16.
      int extra = nWindows % numProc;
      int paddednWindows = nWindows;
      if (extra != 0) {
        paddednWindows = nWindows - extra + numProc;
      }
      numWorkItems = (int) (paddednWindows / numProc);

      if (numProc > 1) {
        logger.info(format(" Number of MPI Processes:  %d", numProc));
        logger.info(format(" Rank of this MPI Process: %d", rank));
        logger.info(format(" Work per process per row: %d", numWorkItems));
      }
    } else {
      world = null;
      numProc = 1;
      rank = 0;
    }

    // Parallelized determination of number of models in BAR calculation.
    maxModelsPJ = new int[numProc][1] // Note: int does not work with PJ gather commands --> int[][]
    buffersMax = new IntegerBuf[numProc];
    for(int p = 0; p < numProc; p++){
      Arrays.fill(maxModelsPJ[p], MIN_VALUE)
      buffersMax[p] = IntegerBuf.buffer(maxModelsPJ[p])
    }

    myBufferMax = buffersMax[rank]
    if (useTinkerBAR) {
      for (int w = 0; w < nWindows; w++) {
        int windowRank = w % numProc;
        if (windowRank == rank) {
          for (String fileName : fullFilePaths[w]) {
            try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {
              String data;
              int xyzCount = 0;
              while ((data = br.readLine()) != null) {
                if (data.contains(".xyz") || data.contains(".pdb")) {
                  xyzCount++
                  String[] tokens = data.trim().split(" +");
                  int numModels = parseInt(tokens[0])
                  if (numModels > maxModelsPJ[rank][0]) {
                    maxModelsPJ[rank][0] = numModels
                  }
                }
                if (xyzCount == numTopologies) {
                  break;
                }
              }
            } catch (IOException fileNotFoundException) {
              logger.warning(format(" Exception reading %s:\n %s", fileName, fileNotFoundException));
            }
          }
        }
      }
    } else {
      // Open archive files for energy evaluation.
      for (int t = 0; t < numTopologies; t++) {
        // Might need parameter/patch information from original file, therefore load original file.
        MolecularAssembly unused = alchemical.openFile(algorithmFunctions, topology, threadsPer, filenames[t], t)
        // Extract systemFilter that will be used to read archives.
        SystemFilter systemFilter = algorithmFunctions.getFilter();
        for (int w = 0; w < nWindows; w++) {
          int windowRank = w % numProc;
          if (windowRank == rank) {
            String fileName = fullFilePaths[w][t];
            // Update assembly to use archive file.
            systemFilter.getActiveMolecularSystem().setFile(new File(fileName));
            // Determine minimum number of assemblies.
            int numModels = systemFilter.countNumModels()
            if (numModels > maxModelsPJ[rank][0]) {
              maxModelsPJ[rank][0] = numModels;
            }
          }
        }
      }
    }
    int maxModels;
    if(useMPI) {
      if(rank == 0 && numProc > 1){
        logger.info(" MPI is only implemented for energy calculation (useTinkerBar or from files).")
      }
      world.allGather(myBufferMax, buffersMax)
      world.barrier()
      int maximum = MIN_VALUE;
      for (int[] numModel : maxModelsPJ) {
        for(int value: numModel){
          if (value > maximum) {
            maximum = value;
          }
        }
      }
      for(int p = 0; p < numProc; p++){
        Arrays.fill(maxModelsPJ[p], maximum)
      }
      maxModels = maxModelsPJ[rank][0]
    }else {
      int maximum = MIN_VALUE;
      for (int[] numModel : maxModelsPJ) {
        for(int value: numModel){
          if (value > maximum) {
            maximum = value;
          }
        }
      }
      maxModels = maximum;
    }
    logger.info(format(" Maximum number of models is: %3d", maxModels));

    // Determine number of snapshots to include in evaluation.
    if(endingSnapshot == 0){
      endingSnapshot = maxModels - 1
    }
    int snapshots = (endingSnapshot - startingSnapshot) + 1;
    if(snapshots > maxModels){
      logger.warning(format(" Specified number of snapshots (%5d) is greater than snapshots available (%5d).", snapshots, maxModels))
    }

    // Initialize objects to be used during BAR.
    energiesLow = new double[nWindows][snapshots]
    energiesAt = new double[nWindows][snapshots]
    energiesHigh = new double[nWindows][snapshots]
    volume = new double[nWindows][snapshots]

    energiesLowPJ = new double[numProc][numWorkItems][snapshots]
    energiesAtPJ = new double[numProc][numWorkItems][snapshots]
    energiesHighPJ = new double[numProc][numWorkItems][snapshots]
    volumePJ = new double[numProc][numWorkItems][snapshots]
    // Default values to impossible so that we know which should be removed in the end.
    //   Otherwise use minimum number of snapshots for all windows.
    for(int i = 0; i < numProc; i++){
      for(int j = 0; j < numWorkItems; j++){
        Arrays.fill(energiesLowPJ[i][j], Double.NaN);
        Arrays.fill(energiesAtPJ[i][j], Double.NaN);
        Arrays.fill(energiesHighPJ[i][j], Double.NaN);
        Arrays.fill(volumePJ[i][j], Double.NaN);
      }
    }

    buffersLow = new DoubleBuf[numProc];
    buffersAt = new DoubleBuf[numProc];
    buffersHigh = new DoubleBuf[numProc];
    buffersVolume = new DoubleBuf[numProc];

    for (int p = 0; p < numProc; p++) {
      buffersLow[p] = DoubleBuf.buffer(energiesLowPJ[p]);
      buffersAt[p] = DoubleBuf.buffer(energiesAtPJ[p])
      buffersHigh[p] = DoubleBuf.buffer(energiesHighPJ[p])
      buffersVolume[p] = DoubleBuf.buffer(volumePJ[p])
    }

    myBufferLow = buffersLow[rank]
    myBufferAt = buffersAt[rank]
    myBufferHigh = buffersHigh[rank]
    myBufferVolume = buffersVolume[rank]

    double[] currentLambdas
    double[][] energy

    int nCurrLambdas
    for (int w = 0; w < nWindows; w++) {
      if (w == 0) {
        currentLambdas = new double[2]
        currentLambdas[0] = lambdaValues[w]
        currentLambdas[1] = lambdaValues[w + 1]
      } else if (w == nWindows - 1) {
        currentLambdas = new double[2]
        currentLambdas[0] = lambdaValues[w - 1]
        currentLambdas[1] = lambdaValues[w]
      } else {
        currentLambdas = new double[3]
        currentLambdas[0] = lambdaValues[w - 1]
        currentLambdas[1] = lambdaValues[w]
        currentLambdas[2] = lambdaValues[w + 1]
      }
      nCurrLambdas = currentLambdas.length;
      energy = new double[nCurrLambdas][]

      if (useTinkerBAR) {
        File barFile = new File(fullFilePaths[w][0])
        barOpeners[w] = new BARFilter(barFile, startingSnapshot, endingSnapshot)
        if (w != nWindows - 1) {
          barOpeners[w].readFile()
        }
        int windowRank = w % numProc;
        if (windowRank == rank) {
          int workIndex = (int) (w / numProc);
          if (w == 0) {
            energiesAtPJ[rank][workIndex] = barOpeners[w].getE1l1()
            energiesHighPJ[rank][workIndex] = barOpeners[w].getE1l2()
          } else if (w == nWindows - 1) {
            energiesLowPJ[rank][workIndex] = barOpeners[w - 1].getE2l1()
            energiesAtPJ[rank][workIndex] = barOpeners[w - 1].getE2l2()
          } else if (w > 0 && w < nWindows - 1) {
            energiesLowPJ[rank][workIndex] = barOpeners[w - 1].getE2l1()
            energiesAtPJ[rank][workIndex] = barOpeners[w].getE1l1()
            energiesHighPJ[rank][workIndex] = barOpeners[w].getE1l2()
          }

          if (isPBC) {
            volumePJ[rank][workIndex] = barOpeners[w].getVolume1()
          }
        }
      } else {
        int windowRank = w % numProc;
        if (windowRank == rank) {
          int workIndex = (int) (w / numProc);
          volumePJ[rank][workIndex] = getEnergyForLambdas(molecularAssemblies, currentLambdas,
                  fullFilePaths[w], energy, isPBC, nSymm)
          if (w == 0) {
            energiesAtPJ[rank][workIndex] = energy[0]
            energiesHighPJ[rank][workIndex] = energy[1]
          } else if (w == nWindows - 1) {
            energiesLowPJ[rank][workIndex] = energy[0]
            energiesAtPJ[rank][workIndex] = energy[1]
          } else if (w > 0 && w < nWindows - 1) {
            energiesLowPJ[rank][workIndex] = energy[0]
            energiesAtPJ[rank][workIndex] = energy[1]
            energiesHighPJ[rank][workIndex] = energy[2]
          }
        }
      }
    }
    // Gather all values from processes and place in object with desired dimensions.
    gatherAllValues();

    // Create file objects to write out TINKER style bar files.
    String tinkerFilePath = ""
    if (tinkerBAR) {
      String tinkerDirectoryPath = directoryPath + File.separator + "barFiles"
      File directory = new File(tinkerDirectoryPath)
      tinkerFilePath = tinkerDirectoryPath + File.separator
      directory.mkdir()
    }

    File xyzFile = new File(filenames.get(0))
    double[] energyMean = new double[nWindows]
    double[] energySD = new double[nWindows]
    double[] energyVar = new double[nWindows]
    // Intensive calculations are done. Utilizing node 0 only for remaining calculations/output.
    if(!useMPI || rank == 0) {
      for (int w = 0; w < nWindows + 1; w++) {
        if (w < nWindows) {
          if (tinkerBAR) {
            if (w == 0) {
              barWriters[w] = new BARFilter(xyzFile, energiesAt[w], energiesHigh[w], energiesLow[w + 1],
                      energiesAt[w + 1], volume[w], volume[w + 1], temperature)
            } else if (w != nWindows - 1) {
              barWriters[w] = new BARFilter(xyzFile, energiesAt[w], energiesHigh[w], energiesLow[w + 1],
                      energiesAt[w + 1], volume[w], volume[w + 1], temperature)
            }
            if (w != nWindows - 1) {
              String barFileName = tinkerFilePath + "energy_" + w.toString() + ".bar"
              barWriters[w].writeFile(barFileName, isPBC)
            }
          }

          BootStrapStatistics energyStats = new BootStrapStatistics(energiesAt[w])
          energyMean[w] = energyStats.mean
          energySD[w] = energyStats.sd
          energyVar[w] = energyStats.var
        }
      }
      // Finish writing BAR files before evaluating windows in case of error evaluating a window.
      for (int w = 0; w < nWindows + 1; w++) {
        if (w == nWindows) {
          logger.info("\n\n Evaluating Overall:")
        } else {
          logger.info(format("\n\n Evaluating Window %d:", w))
        }

        if (w == 0) {
          currentLambdas = new double[2]
          currentLambdas[0] = lambdaValues[w]
          currentLambdas[1] = lambdaValues[w + 1]
        } else if (w == nWindows - 1) {
          currentLambdas = new double[2]
          currentLambdas[0] = lambdaValues[w - 1]
          currentLambdas[1] = lambdaValues[w]
        } else if (w == nWindows) {
          currentLambdas = lambdaValues
        } else {
          currentLambdas = new double[3]
          currentLambdas[0] = lambdaValues[w - 1]
          currentLambdas[1] = lambdaValues[w]
          currentLambdas[2] = lambdaValues[w + 1]
        }

        nCurrLambdas = currentLambdas.length;
        double[][] energyWindowLow = new double[nCurrLambdas][]
        double[][] energyWindowAt = new double[nCurrLambdas][]
        double[][] energyWindowHigh = new double[nCurrLambdas][]

        if (w == 0) {
          energyWindowLow[0] = energiesLow[w]
          energyWindowLow[1] = energiesLow[w + 1]

          energyWindowAt[0] = energiesAt[w]
          energyWindowAt[1] = energiesAt[w + 1]

          energyWindowHigh[0] = energiesHigh[w]
          energyWindowHigh[1] = energiesHigh[w + 1]
        } else if (w == nWindows - 1) {
          energyWindowLow[0] = energiesLow[w - 1]
          energyWindowLow[1] = energiesLow[w]

          energyWindowAt[0] = energiesAt[w - 1]
          energyWindowAt[1] = energiesAt[w]

          energyWindowHigh[0] = energiesHigh[w - 1]
          energyWindowHigh[1] = energiesHigh[w]
        } else if (w == nWindows) {
          energyWindowLow = energiesLow
          energyWindowAt = energiesAt
          energyWindowHigh = energiesHigh
        } else {
          energyWindowLow[0] = energiesLow[w - 1]
          energyWindowLow[1] = energiesLow[w]
          energyWindowLow[2] = energiesLow[w + 1]

          energyWindowAt[0] = energiesAt[w - 1]
          energyWindowAt[1] = energiesAt[w]
          energyWindowAt[2] = energiesAt[w + 1]

          energyWindowHigh[0] = energiesHigh[w - 1]
          energyWindowHigh[1] = energiesHigh[w]
          energyWindowHigh[2] = energiesHigh[w + 1]
        }

        SequentialEstimator bar = new BennettAcceptanceRatio(currentLambdas, energyWindowLow,
                energyWindowAt, energyWindowHigh, new double[]{temperature}, eps, nIterations)
        SequentialEstimator forwards = bar.getInitialForwardsGuess()
        SequentialEstimator backwards = bar.getInitialBackwardsGuess()

        EstimateBootstrapper barBS = new EstimateBootstrapper(bar)
        EstimateBootstrapper forBS = new EstimateBootstrapper(forwards)
        EstimateBootstrapper backBS = new EstimateBootstrapper(backwards)

        int volumeLength = volume.length;
        long bootstrap = min(MAX_BOOTSTRAP_TRIALS, min(volumeLength, volumeLength))
        if (w == nWindows) {
          logger.info("\n Free Energy Difference:\n")
        } else {
          logger.info(format("\n Free Energy Difference for Window %d\n", w))
        }

        long time = -System.nanoTime()
        forBS.bootstrap(bootstrap)
        time += System.nanoTime()
        logger.fine(format(" Forward FEP Bootstrap Complete:      %7.4f sec", time * NS2SEC))
        forwardFEP = forBS.getTotalFE()
        forwardEnthalpy = forBS.getTotalEnthalpy()
        double varForeFE = forBS.getTotalUncertainty()
        double varEnthalpyFore = forBS.getTotalEnthalpyUncertainty()
        logger.info(format(" Free energy via Forwards FEP:   %12.4f +/- %6.4f kcal/mol.", forwardFEP, varForeFE))

        time = -System.nanoTime()
        backBS.bootstrap(bootstrap)
        time += System.nanoTime()
        logger.fine(format(" Backward FEP Bootstrap Complete:     %7.4f sec", time * NS2SEC))

        backwardFEP = backBS.getTotalFE()
        backwardEnthalpy = backBS.getTotalEnthalpy()
        double varBackFE = backBS.getTotalUncertainty()
        double varEnthalpyBack = backBS.getTotalEnthalpyUncertainty()
        logger.info(format(" Free energy via Backwards FEP:  %12.4f +/- %6.4f kcal/mol.", backwardFEP, varBackFE))
        barEnergy = bar.getFreeEnergy()

        logger.info(format(" Free energy via BAR Iteration:  %12.4f +/- %6.4f kcal/mol.", barEnergy, bar.getUncertainty()))
        time = -System.nanoTime()
        barBS.bootstrap(bootstrap)
        time += System.nanoTime()
        logger.fine(format(" BAR Bootstrap Complete:              %7.4f sec", time * NS2SEC))

        barEnergyBS = barBS.getTotalFE()
        double varBARFE = barBS.getTotalUncertainty()
        barEnthalpy = barBS.getTotalEnthalpy()
        double varEnthalpy = barBS.getTotalEnthalpyUncertainty()
        logger.info(format(" Free energy via BAR Bootstrap:  %12.4f +/- %6.4f kcal/mol.", barEnergyBS, varBARFE))

        if (w == nWindows) {
          logger.info("\n Enthalpy from Potential Energy Averages:\n")
          for (int n = 0; n < nWindows; n++) {
            logger.info(format(" Average Energy for State %d:       %12.4f +/- %6.4f kcal/mol.", n, energyMean[n], energySD[n]))
          }
          double enthalpyDiff = energyMean[nWindows - 1] - energyMean[0]
          double enthalpyDiffSD = Math.sqrt(energyVar[nWindows - 1] + energyVar[0])
          logger.info(format(" Enthalpy via Direct Estimate:     %12.4f +/- %6.4f kcal/mol.", enthalpyDiff, enthalpyDiffSD))
          logger.info("\n Enthalpy and Entropy:\n")
        } else {
          logger.info(format("\n Enthalpy and Entropy for Window %d\n", w))
        }

        forwardEntropy = (forwardEnthalpy - forwardFEP) / temperature
        backwardEntropy = (backwardEnthalpy - backwardFEP) / temperature

        logger.info(format(" Enthalpy via Forward FEP:       %12.4f +/- %6.4f kcal/mol.", forwardEnthalpy, varEnthalpyFore))
        logger.info(format(" Entropy via Forward FEP:        %12.4f kcal/mol/K.", forwardEntropy))
        logger.info(format(" Forward FEP -T*ds Value:        %12.4f kcal/mol.", -(forwardEntropy * temperature)))

        logger.info(format("\n Enthalpy via Backward FEP:      %12.4f +/- %6.4f kcal/mol.", backwardEnthalpy, varEnthalpyBack))
        logger.info(format(" Entropy via Backward FEP:       %12.4f kcal/mol/K.", backwardEntropy))
        logger.info(format(" Backward FEP -T*ds Value:       %12.4f kcal/mol.", -(backwardEntropy * temperature)))

        double tsBar = barEnthalpy - barEnergyBS
        double sBAR = tsBar / (temperature)
        logger.info(format("\n Enthalpy via BAR:               %12.4f +/- %6.4f kcal/mol.", barEnthalpy, varEnthalpy))
        logger.info(format(" Entropy via BAR:                %12.4f kcal/mol/K.", sBAR))
        logger.info(format(" BAR Estimate of -T*ds:          %12.4f kcal/mol.", -(tsBar)))
      }
    }

    return this

  }

  /**
   * This method calls <code>world.gather</code> to collect numProc values.
   */
  private void gatherAllValues() {
    if (useMPI) {
      try {
        world.gather(0, myBufferLow, buffersLow);
        world.gather(0, myBufferAt, buffersAt);
        world.gather(0, myBufferHigh, buffersHigh);
        world.gather(0, myBufferVolume, buffersVolume);
        if (rank == 0) {
          for (int workItem = 0; workItem < numWorkItems; workItem++) {
            for (int proc = 0; proc < numProc; proc++) {
              final int index = numProc * workItem + proc;
              // Do not include padded results.
              if (index < nWindows) {
                energiesLow[index] = energiesLowPJ[proc][workItem];
                energiesAt[index] = energiesAtPJ[proc][workItem];
                energiesHigh[index] = energiesHighPJ[proc][workItem];
                volume[index] = volumePJ[proc][workItem];
              }
            }
          }
        } else {
          for (int workItem = 0; workItem < numWorkItems; workItem++) {
            final int index = numProc * workItem + rank;
            // Do not include padded results.
            if (index < nWindows) {
              energiesLow[index] = energiesLowPJ[rank][workItem];
              energiesAt[index] = energiesAtPJ[rank][workItem];
              energiesHigh[index] = energiesHighPJ[rank][workItem];
              volume[index] = volumePJ[rank][workItem];
            }
          }
        }
      } catch (Exception ex) {
        logger.severe(" Exception collecting distance values." + ex + Utilities.stackTraceToString(ex));
      }
    } else {
      for (int i = 0; i < nWindows; i++) {
        energiesLow[i] = energiesLowPJ[rank][i];
        energiesAt[i] = energiesAtPJ[rank][i];
        energiesHigh[i] = energiesHighPJ[rank][i];
        volume[i] = volumePJ[rank][i];
      }
    }
  }

  private double[] getEnergyForLambdas(MolecularAssembly[] topologies, double[] lambdaValues,
                                       String[] arcFileName, double[][] energy, boolean isPBC, int nSymm) {
    for (int j = 0; j < numTopologies; j++) {
      File archiveFile = new File(arcFileName[j])
      openers[j].setFile(archiveFile)
      topologies[j].setFile(archiveFile)
      StringBuilder sb = new StringBuilder(format(
              "\n Evaluating energies for %s\n ", arcFileName[j]))
      logger.info(sb as String)
    }
    int nSnapshots = openers[0].countNumModels()

    double[] x = new double[potential.getNumberOfVariables()]
    double[] vol = new double[nSnapshots]
    int nLambdas = lambdaValues.length
    for (int k = 0; k < nLambdas; k++) {
      energy[k] = new double[nSnapshots]
    }

    LambdaInterface linter1 = (LambdaInterface) potential

    int endWindow = nWindows - 1
    String endWindows = endWindow + File.separator

    if (arcFileName[0].contains(endWindows)) {
      logger.info(format(" %s     %s   %s     %s   %s ", "Snapshot", "Lambda Low",
              "Energy Low", "Lambda At", "Energy At"))
    } else if (arcFileName[0].contains("0/")) {
      logger.info(format(" %s     %s   %s     %s   %s ", "Snapshot", "Lambda At",
              "Energy At", "Lambda High", "Energy High"))
    } else {
      logger.info(format(" %s     %s   %s     %s   %s     %s   %s ", "Snapshot", "Lambda Low",
              "Energy Low", "Lambda At", "Energy At", "Lambda High", "Energy High"))
    }

    for (int i = 0; i < nSnapshots; i++) {
      boolean resetPosition = (i == 0)
      int nOpeners = openers.length
      for (int n = 0; n < nOpeners; n++) {
        openers[n].readNext(resetPosition, false)
      }

      x = potential.getCoordinates(x)
      nLambdas = lambdaValues.length
      for (int k = 0; k < nLambdas; k++) {
        double lambda = lambdaValues[k]
        linter1.setLambda(lambda)
        // Consider falling back to DIRECT if SCF.
        energy[k][i] = potential.energy(x, false)
      }

      if (nLambdas == 2) {
        logger.info(format(" %8d     %6.3f   %14.4f     %6.3f   %14.4f ", i + 1, lambdaValues[0],
                energy[0][i], lambdaValues[1], energy[1][i]))
      } else {
        logger.info(format(" %8d     %6.3f   %14.4f     %6.3f   %14.4f     %6.3f   %14.4f ", i + 1,
                lambdaValues[0],
                energy[0][i], lambdaValues[1], energy[1][i], lambdaValues[2], energy[2][i]))
      }

      if (isPBC) {
        Crystal unitCell = potential.getCrystal().getUnitCell()
        vol[i] = unitCell.volume / nSymm
        logger.info(format(" %8d %14.4f",
                i + 1, vol[i]))
      }
    }

    return vol
  }

/**
 * {@inheritDoc}
 */
  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (potential == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = new ArrayList<>()
      if (potential != null) {
        potentials.add(potential)
      }
    }
    return potentials
  }


  double getBarEnergy() {
    return barEnergy
  }

  double getBarEnergyBS() {
    return barEnergyBS
  }

  double getFepFor() {
    return forwardFEP
  }

  double getFepBack() {
    return backwardFEP
  }

  double gethFor() {
    return forwardEnthalpy
  }

  double gethBack() {
    return backwardEnthalpy
  }

  double gethBAR() {
    return barEnthalpy
  }

  double getsFor() {
    return forwardEntropy
  }

  double getsBack() {
    return backwardEntropy
  }

  double getsBAR() {
    return barEntropy
  }

}
