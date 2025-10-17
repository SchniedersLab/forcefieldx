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
package ffx.algorithms.commands;

import ffx.algorithms.ParallelStateEnergy;
import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.numerics.Potential;
import ffx.numerics.estimator.FreeEnergyDifferenceReporter;
import ffx.potential.MolecularAssembly;
import ffx.potential.cli.AlchemicalOptions;
import ffx.potential.cli.TopologyOptions;
import ffx.potential.parsers.BARFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.utilities.FFXBinding;
import org.apache.commons.configuration2.Configuration;
import org.apache.commons.io.FilenameUtils;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import static java.lang.String.format;

/**
 * The BAR script finds the free energy difference across a lambda window. It presently assumes
 * that the number of files composing the first end of the window equals the number of files
 * composing the other end.
 * <br>
 * Usage:
 * <br>
 * ffxc BAR [options] &lt;structures1&gt; &lt;structures2&gt;
 */
@Command(description = " Evaluates a free energy change with the Bennett Acceptance Ratio algorithm using pregenerated snapshots.", name = "BAR")
public class BAR extends AlgorithmsCommand {

  @Mixin
  private AlchemicalOptions alchemicalOptions;

  @Mixin
  private TopologyOptions topologyOptions;

  /**
   * --eb or --evaluateBARFiles Read in Tinker BAR files and evaluate the free energy difference.
   */
  @Option(names = {"--eb", "--evaluateBAR"}, paramLabel = "false", defaultValue = "false",
      description = "Evaluate Free Energy Differences from Tinker BAR files.")
  private boolean evaluateBARFiles = false;

  @Option(names = {"--l2", "--lambdaTwo"}, paramLabel = "1.0",
      description = "Lambda value for the upper edge of the window")
  private double lambda2 = 1.0;

  @Option(names = {"-t", "--temperature"}, paramLabel = "298.15",
      description = "Temperature for system")
  private double temperature = 298.15;

  @Option(names = {"--dV", "--volume"}, paramLabel = "false",
      description = "Write out snapshot volumes to the Tinker BAR file.")
  private boolean includeVolume = false;

  @Option(names = {"--ns", "--nStates"}, paramLabel = "2",
      description = "If not equal to two, auto-determine lambda values and subdirectories (overrides other flags).")
  private int nStates = 2;

  @Option(names = {"--sa", "--sortedArc"}, paramLabel = "false",
      description = "If set, use sorted archive values.")
  private boolean sortedArc = false;

  @Option(names = {"--ss", "--startSnapshot"}, paramLabel = "0",
      description = "Start at this snapshot when reading in Tinker BAR files (indexed from 0).")
  private int startingSnapshot = 0;

  @Option(names = {"--es", "--endSnapshot"}, paramLabel = "0",
      description = "End at this snapshot when reading in Tinker BAR files (indexed from 0).")
  private int endingSnapshot = 0;

  /**
   * --ni or --nIterations Maximum number of allowable iterations for BAR calculation.
   */
  @Option(names = {"--ni", "--nIterations"}, paramLabel = "1000",
      description = "Specify the maximum number of iterations for BAR convergence.")
  private int nIterations = 1000;

  /**
   * -e or --eps Convergence criterion for BAR iteration.
   */
  @Option(names = {"-e", "--eps"}, paramLabel = "1.0E-4",
      description = "Specify convergence cutoff for BAR calculation.")
  private double eps = 1.0E-4;

  /**
   * The final argument(s) can be filenames for lambda windows in order.
   */
  @Parameters(arity = "0..*", paramLabel = "files",
      description = "A single PDB/XYZ when windows are auto-determined (or two for dual topology). Two trajectory files for BAR between two ensembles (or four for dual topology).")
  private List<String> filenames = null;

  /**
   * Additional properties.
   */
  private Configuration additionalProperties;

  /**
   * Molecular assemblies to compute energies.
   */
  MolecularAssembly[] topologies;

  /**
   * Number of Topologies.
   */
  private int numTopologies;

  /**
   * Potential object for the crystal.
   */
  CrystalPotential potential;

  /**
   * Number of input files.
   */
  int nFiles;

  /**
   * Input files for BAR.
   */
  private String[] files;

  /**
   * Lambda values for each state.
   */
  double[] lambdaValues;

  /**
   * Compute and log out the free energy differences.
   */
  private FreeEnergyDifferenceReporter reporter;

  /**
   * BAR Constructor.
   */
  public BAR() {
    super();
  }

  /**
   * BAR Constructor.
   *
   * @param binding The Binding to use.
   */
  public BAR(FFXBinding binding) {
    super(binding);
  }

  /**
   * BAR constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public BAR(String[] args) {
    super(args);
  }

  /**
   * Sets an optional Configuration with additional properties.
   *
   * @param additionalProps Additional properties configuration.
   */
  public void setProperties(Configuration additionalProps) {
    this.additionalProperties = additionalProps;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BAR run() {
    // Initialize the script.
    if (!init()) {
      return this;
    }

    /**
     * Load user supplied files into an array.
     */
    if (filenames == null) {
      nFiles = 0;
      files = null;
    } else {
      nFiles = filenames.size();
      files = new String[nFiles];
      for (int i = 0; i < nFiles; i++) {
        files[i] = filenames.get(i);
      }
    }

    // Evaluate free energy differences from BAR files.
    if (evaluateBARFiles) {
      // 1) Define the number of states and concordant BAR files (windows = states - 1).
      // 2) Read in the energy and volume data from BAR files.
      // 3) Compute the free energy differences.
      if (nFiles == 0) {
        logger.info(" Please supply a BAR file or directory to be evaluated.");
        return this;
      }

      BARFilter[] barFilters = null;
      if (nStates == 2 && nFiles == 1) {
        // If the number of states is 2, evaluate a single BAR file.
        logger.info(format(" Evaluating a single BAR file: %s.", files[0]));
        barFilters = new BARFilter[1];
        barFilters[0] = new BARFilter(new File(files[0]), startingSnapshot, endingSnapshot);
        // Define the lambda values.
        lambdaValues = new double[2];
        lambdaValues[0] = alchemicalOptions.getInitialLambda();
        lambdaValues[1] = lambda2;
      } else if (nStates > 2 && nFiles == 1) {
        // If the number of states is greater than 2, evaluate multiple BAR files from a directory.
        logger.info(format(" Evaluating BAR files from directory %s.", files[0]));
        // Auto-determine bar files from the "bar" subdirectory.
        logger.info(" Auto-detecting BAR files and lambda values.");
        Path subdirectoryPath = Paths.get(files[0]);
        // Regular expression to match files with the form "filename_XX.bar"
        Pattern pattern = Pattern.compile(".*_(\\d+)\\.bar");
        // List to hold file paths
        List<String> fileList = new ArrayList<>();
        try (Stream<Path> paths = Files.list(subdirectoryPath)) {
          paths.filter(Files::isRegularFile)
              .filter(path -> {
                Matcher matcher = pattern.matcher(path.getFileName().toString());
                return matcher.matches();
              })
              .sorted(Comparator.comparingInt(path -> {
                Matcher matcher = pattern.matcher(path.getFileName().toString());
                matcher.matches();
                return Integer.parseInt(matcher.group(1));
              }))
              .forEach(path -> fileList.add(path.toString()));
        } catch (IOException e) {
          logger.info(" Error reading files from a directory.\n" + e.toString());
          e.printStackTrace();
          return this;
        }

        int numFiles = fileList.size();
        if (numFiles != nStates - 1) {
          logger.info(format(" The number of bar files (%d) does not concord with the number of states (%d).", numFiles, nStates));
          return this;
        }

        // Define the lambda values.
        lambdaValues = new double[nStates];
        for (int l = 0; l < nStates; l++) {
          lambdaValues[l] = alchemicalOptions.getInitialLambda(nStates, l, true);
        }
        barFilters = new BARFilter[numFiles];
        for (int i = 0; i < numFiles; i++) {
          String filename = fileList.get(i);
          File barFile = new File(filename);
          barFilters[i] = new BARFilter(barFile, startingSnapshot, endingSnapshot);
          logger.info(format(" Window L=%6.4f to L=%6.4f from %s",
              lambdaValues[i], lambdaValues[i + 1], filename));
        }
      }

      // Allocate space for energy and volume arrays.
      double[][] energiesLow = new double[nStates][];
      double[][] energiesAt = new double[nStates][];
      double[][] energiesHigh = new double[nStates][];
      double[][] volume = new double[nStates][];
      double[] temperatures = new double[nStates];

      // Load BAR files.
      readBARFiles(barFilters, energiesLow, energiesAt, energiesHigh, volume, temperatures);

      // Compute free energy differences.
      reporter = new FreeEnergyDifferenceReporter(nStates, lambdaValues, temperatures, eps, nIterations,
          energiesLow, energiesAt, energiesHigh, volume);
      reporter.report();

      // Exit the script.
      return this;
    }

    // Create BAR files for later evaluation of free energy differences.
    logger.info(" Writing BAR files.");

    numTopologies = 1;
    if (nFiles <= 1 && nStates < 2) {
      logger.info(" At least two states must be specified");
      return this;
    } else if (nFiles == 1 && nStates >= 2) {
      logger.info(format(" Auto-detecting %d states for single topology:\n %s.", nStates, files[0]));
    } else if (nFiles == 2 && nStates >= 2) {
      logger.info(format(" Auto-detecting %d states for dual topology:\n %s\n %s.", nStates, files[0], files[1]));
      numTopologies = 2;
    } else if (nFiles == 2) {
      logger.info(format(" Applying BAR between two single topology ensembles:\n %s\n %s.", files[0], files[1]));
      nStates = 2;
    } else if (nFiles == 4) {
      logger.info(format(" Applying BAR between two dual topology ensembles:\n %s %s\n %s %s.",
          files[0], files[1], files[2], files[3]));
      numTopologies = 2;
      nStates = 2;
    } else {
      logger.info(format(" Inconsistent input of files (%3d) and/or states (%3d).", nFiles, nStates));
      return this;
    }

    boolean autodetect = false;
    if (nStates > 2) {
      autodetect = true;
      lambdaValues = new double[nStates];
      for (int i = 0; i < nStates; i++) {
        lambdaValues[i] = alchemicalOptions.getInitialLambda(nStates, i, true);
      }
    } else {
      // Otherwise we assume two ensembles at the given lambda values.
      lambdaValues = new double[2];
      lambdaValues[0] = alchemicalOptions.getInitialLambda(2, 0, true);
      lambdaValues[1] = lambda2;
      nStates = 2;
    }

    // Could set "getInitialLambda"'s quiet flag to false, but better logging here?
    logger.info(" Lambda values for each window: ");
    int nLambda = lambdaValues.length;
    for (int i = 0; i < nLambda; i++) {
      logger.info(format(" Window %3d: %6.4f", i, lambdaValues[i]));
    }

    // Open assemblies if not using Tinker BAR files.
    boolean isPBC;

    // Allocate space for each topology
    int threadsPerTopology = topologyOptions.getThreadsPerTopology(numTopologies);
    topologies = new MolecularAssembly[numTopologies];
    SystemFilter[] openers = new SystemFilter[numTopologies];

    alchemicalOptions.setAlchemicalProperties();
    topologyOptions.setAlchemicalProperties(numTopologies);

    if (numTopologies == 2) {
      logger.info(format(" Initializing two topologies for each window."));
    } else {
      logger.info(format(" Initializing a single topology for each window."));
    }

    for (int i = 0; i < numTopologies; i++) {
      MolecularAssembly ma = alchemicalOptions.openFile(algorithmFunctions, topologyOptions,
          threadsPerTopology, filenames.get(i), i);
      topologies[i] = ma;
      openers[i] = algorithmFunctions.getFilter();
    }

    StringBuilder sb = new StringBuilder(format(
        "\n Using BAR to analyze a free energy change for %s\n ", filenames));
    potential = (CrystalPotential) topologyOptions.assemblePotential(topologies, sb);
    Crystal unitCell = potential.getCrystal().getUnitCell();
    isPBC = includeVolume && !unitCell.aperiodic();

    String[][] fullFilePaths;
    String directoryPath;
    if (nFiles > 0) {
      // For Dual-Topology systems that only compare two states, simplify filePaths to use only the 2 files belonging
      // to each ensemble
      int dtIndex = 0;
      if (numTopologies == 2 && nFiles == 4) {
        fullFilePaths = new String[nStates][2];
        dtIndex = 2;
      } else {
        fullFilePaths = new String[nStates][nFiles];
      }

      // Loop over ensembles
      for (int i = 0; i < nStates; i++) {
        // Loop over user supplied files.
        for (int j = 0; j < nFiles; j++) {
          // For Dual-Topology, start at index 2 for second ensemble
          if (i == 1 && j == 0) {
            j = dtIndex;
          }
          File file = new File(files[j]);
          directoryPath = file.getAbsoluteFile().getParent() + File.separator;
          String archiveName;
          if (sortedArc) {
            archiveName = FilenameUtils.getBaseName(files[j]) + "_E" + String.valueOf(i) + ".arc";
          } else {
            archiveName = FilenameUtils.getBaseName(files[j]) + ".arc";
          }
          if (!autodetect) {
            // Path to a file in the same directory as supplied archives.
            fullFilePaths[i][j - dtIndex] = directoryPath + File.separator + archiveName;
          } else {
            // Paths to auto-detected subdirectories.
            fullFilePaths[i][j - dtIndex] = directoryPath + i + File.separator + archiveName;
          }
          // For Dual-Topology, stop after two files for first ensemble
          if (i == 0 && j == 1) {
            j += dtIndex * nFiles;
          }
        }
      }

      // Initialize the ParallelEnergy class.
      ParallelStateEnergy pe = new ParallelStateEnergy(nStates, lambdaValues,
          topologies, potential, fullFilePaths, openers);

      // Evaluate energies from snapshots.
      double[][] energiesLow = new double[nStates][];
      double[][] energiesAt = new double[nStates][];
      double[][] energiesHigh = new double[nStates][];
      double[][] volume = new double[nStates][];
      pe.evaluateStates(energiesLow, energiesAt, energiesHigh, volume);

      // Only node 0 will write out Tinker BAR files.
      if (pe.getRank() != 0) {
        return this;
      }

      // Create file objects to write out TINKER style bar files.
      File file = new File(files[0]);
      directoryPath = file.getAbsoluteFile().getParent() + File.separator;
      String tinkerDirectoryPath = directoryPath + File.separator + "windows";
      File directory = new File(tinkerDirectoryPath);
      String barFilePath = tinkerDirectoryPath + File.separator;
      directory.mkdir();

      for (int state = 0; state < nStates - 1; state++) {
        File xyzFile = new File(filenames.get(0));
        BARFilter barFilter;
        barFilter = new BARFilter(xyzFile, energiesAt[state], energiesHigh[state], energiesLow[state + 1],
            energiesAt[state + 1], volume[state], volume[state + 1], temperature, temperature);
        String barFileName = barFilePath + "window_" + String.valueOf(state) + ".bar";
        // Write out the TINKER style bar file. An existing file will be overwritten.
        barFilter.writeFile(barFileName, isPBC, false);
      }
    }

    return this;
  }

  /**
   * Read energy and volume values from BAR files. The first dimension of the energy arrays
   * corresponds to the state index. The second dimension corresponds to the snapshot index.
   * <p>
   * The first dimension of each array should be allocated to the number of states prior
   * to calling this method.
   *
   * @param barFilters   The BAR filters to use.
   * @param energyLow    The energy of each snapshot in state L evaluated at L-dL.
   * @param energyAt     The energy of each snapshot in L evaluated at L.
   * @param energyHigh   The energy from state L evaluated at L+dL.
   * @param volume       The volume of each snapshot from state L.
   * @param temperatures The temperatures for each state.
   */
  void readBARFiles(BARFilter[] barFilters,
                    double[][] energyLow, double[][] energyAt,
                    double[][] energyHigh, double[][] volume,
                    double[] temperatures) {

    // Read energy values from BAR files.
    for (BARFilter barFilter : barFilters) {
      barFilter.readFile();
    }

    // Assign values to the state arrays.
    for (int state = 0; state < nStates; state++) {
      if (state == 0) {
        energyAt[state] = barFilters[state].getE1l1();
        energyHigh[state] = barFilters[state].getE1l2();
        volume[state] = barFilters[state].getVolume1();
        temperatures[state] = barFilters[state].getTemperature1();
      } else if (state == nStates - 1) {
        energyLow[state] = barFilters[state - 1].getE2l1();
        energyAt[state] = barFilters[state - 1].getE2l2();
        volume[state] = barFilters[state - 1].getVolume2();
        temperatures[state] = barFilters[state - 1].getTemperature2();
      } else {
        energyLow[state] = barFilters[state - 1].getE2l1();
        energyAt[state] = barFilters[state].getE1l1();
        energyHigh[state] = barFilters[state].getE1l2();
        volume[state] = barFilters[state].getVolume1();
        temperatures[state] = barFilters[state].getTemperature1();
      }
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public List<Potential> getPotentials() {
    List<Potential> potentials;
    if (potential == null) {
      potentials = Collections.emptyList();
    } else {
      potentials = new ArrayList<>();
      potentials.add(potential);
    }
    return potentials;
  }

  /**
   * Obtain the Free Energy Difference reporter for this class.
   *
   * @return The Free Energy Difference reporter.
   */
  public FreeEnergyDifferenceReporter getReporter() {
    return reporter;
  }

}
