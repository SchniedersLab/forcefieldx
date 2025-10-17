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
package ffx.algorithms.commands.test;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.cli.BarostatOptions;
import ffx.algorithms.thermodynamics.HistogramData;
import ffx.crystal.CrystalPotential;
import ffx.numerics.estimator.BennettAcceptanceRatio;
import ffx.numerics.estimator.EstimateBootstrapper;
import ffx.numerics.estimator.Zwanzig;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.cli.AlchemicalOptions;
import ffx.potential.cli.TopologyOptions;
import ffx.potential.parsers.SystemFilter;
import ffx.utilities.Constants;
import ffx.utilities.FFXBinding;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.configuration2.Configuration;
import org.apache.commons.io.FilenameUtils;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;
import java.util.logging.Level;

import static java.lang.String.format;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.round;

/**
 * The MostBar script uses a single set of archive file(s) from a Metropolized
 * Orthogonal Space Tempering run to evaluate free energy via the Bennett Acceptance Ratio
 * <br>
 * Usage:
 * <br>
 * ffxc test.MostBar [options] &lt;structures1&gt;
 */
@Command(description = " Evaluates free energy of an M-OST run using the BAR estimator.", name = "test.MostBar")
public class MostBar extends AlgorithmsCommand {

  @Mixin
  private AlchemicalOptions alchemicalOptions;

  @Mixin
  private TopologyOptions topologyOptions;

  @Mixin
  private BarostatOptions barostat;

  /**
   * -t or --temperature sets the temperature in Kelvins.
   */
  @Option(names = {"-t", "--temperature"}, paramLabel = "298.15", defaultValue = "298.15",
      description = "Temperature in Kelvins")
  private double temp;

  /**
   * --his or --histogram manually sets the histogram file to read.
   */
  @Option(names = {"--his", "--histogram"}, paramLabel = "file.his", defaultValue = "",
      description = "Manually provided path to a histogram file (otherwise, attempts to autodetect from same directory as input files).")
  private String histogramName;

  /**
   * --lb or --lambdaBins manually specifies a number of evenly spaced lambda bins rather than reading a histogram file.
   */
  @Option(names = {"--lb", "--lambdaBins"}, paramLabel = "autodetected", defaultValue = "-1",
      description = "Manually specified number of lambda bins (else auto-detected from histogram")
  private int lamBins;

  /**
   * -s or --start sets the first frame to be read (usually of the entire archive: first frame of lambda X if --lambdaSorted is set).
   */
  @Option(names = {"-s", "--start"}, paramLabel = "1", defaultValue = "1",
      description = "First snapshot to evaluate (1-indexed, inclusive).")
  private int startFrame;

  /**
   * -s or --start sets the last frame to be read (usually of the entire archive: last frame of lambda X if --lambdaSorted is set).
   */
  @Option(names = {"--fi", "--final"}, paramLabel = "-1", defaultValue = "-1",
      description = "Last snapshot to evaluate (1-indexed, inclusive); leave negative to analyze to end of trajectory.")
  private int finalFrame;

  /**
   * --st or --stride sets the frequency with which snapshots should be evaluated.
   */
  @Option(names = {"--st", "--stride"}, paramLabel = "1", defaultValue = "1",
      description = "First snapshot to evaluate (1-indexed).")
  private int stride;

  /**
   * --bo or --bootstrap sets the number of bootstrap cycles to run.
   */
  @Option(names = {"--bo", "--bootstrap"}, paramLabel = "AUTO", defaultValue = "-1",
      description = "Use this many bootstrap trials to estimate dG and uncertainty; default is 200-100000 (depending on number of frames).")
  private long bootstrap;

  /**
   * --lambdaSorted indicates that this is not a standard M-OST archive, but rather a concatenation of fixed-lambda sampling (affects -s and --fi).
   */
  @Option(names = {"--lambdaSorted"}, paramLabel = "false", defaultValue = "false",
      description = "Input is sorted by lambda rather than simulation progress (sets -s to skip N-1 frames at each lambda value rather than N-1 of all frames).")
  private boolean lambdaSorted;

  /**
   * -v or --verbose enables extra logging (e.g. energies collected, more frequent bootstrap progress updates, etc).
   */
  @Option(names = {"-v", "--verbose"}, paramLabel = "false", defaultValue = "false",
      description = "Print out extra information (e.g. collection of potential energies).")
  private boolean verbose;

  /**
   * The final argument(s) should be filenames for lambda windows in order.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = "Trajectory files for the first end of the window, followed by trajectories for the other end")
  private List<String> filenames;

  private MolecularAssembly[] topologies;
  private SystemFilter[] openers;
  private CrystalPotential potential;
  private LambdaInterface linter;
  private Configuration additionalProperties;

  private List<List<Double>> energiesL;
  private List<List<Double>> energiesUp;
  private List<List<Double>> energiesDown;

  private double[] lamPoints;
  private int[] observations;
  private double lamSep;
  private double halfLamSep;
  private double[] x;
  private final double[] lastEntries = new double[3];
  private static final String energyFormat = "%11.4f kcal/mol";
  private static final String nanFormat = format("%20s", "N/A");
  // First frame (0-indexed).
  private int start;
  // Last frame (0-indexed, exclusive.
  private int end;
  private Level standardLogging = Level.FINE;

  // Interval between logging which bootstrap cycle it is.
  private static final long BOOTSTRAP_PRINT = 50L;
  // Lower/upper bounds for autodetected bootstrap length.
  private static final long MIN_BOOTSTRAP_TRIALS = 200L;
  private static final long MAX_BOOTSTRAP_TRIALS = 50000L;
  private static final long AUTO_BOOTSTRAP_NUMERATOR = 10000000L;

  public void setProperties(CompositeConfiguration addedProperties) {
    additionalProperties = addedProperties;
  }

  /**
   * MostBar Constructor.
   */
  public MostBar() {
    super();
  }

  /**
   * MostBar Constructor.
   *
   * @param binding The Binding to use.
   */
  public MostBar(FFXBinding binding) {
    super(binding);
  }

  /**
   * MostBar constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public MostBar(String[] args) {
    super(args);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public MostBar run() {
    // Begin boilerplate code.
    if (!init()) {
      return this;
    }

    // Determine the number of topologies to be read and allocate the array.
    int numTopologies = topologyOptions.getNumberOfTopologies(filenames);
    int threadsPerTopology = topologyOptions.getThreadsPerTopology(numTopologies);
    topologies = new MolecularAssembly[numTopologies];
    openers = new SystemFilter[numTopologies];

    // Turn on computation of lambda derivatives if softcore atoms exist.
    alchemicalOptions.setAlchemicalProperties();
    topologyOptions.setAlchemicalProperties(numTopologies);

    standardLogging = verbose ? Level.INFO : Level.FINE;

    logger.info(format(" Initializing %d topologies", numTopologies));

    // Read in files.
    if (filenames == null || filenames.isEmpty()) {
      activeAssembly = getActiveAssembly(null);
      if (activeAssembly == null) {
        logger.info(helpString());
        return this;
      }
      filenames = new ArrayList<>();
      filenames.add(activeAssembly.getFile().getName());
      topologies[0] = alchemicalOptions.processFile(topologyOptions, activeAssembly, 0);
    } else {
      logger.info(format(" Initializing %d topologies...", numTopologies));
      for (int i = 0; i < numTopologies; i++) {
        topologies[i] = alchemicalOptions.openFile(algorithmFunctions,
            topologyOptions, threadsPerTopology, filenames.get(i), i);
        openers[i] = algorithmFunctions.getFilter();
      }
    }

    StringBuilder sb = new StringBuilder("\n Using BAR to analyze an M-OST free energy change for systems ");
    potential = (CrystalPotential) topologyOptions.assemblePotential(topologies, sb);
    potential = barostat.checkNPT(topologies[0], potential);
    linter = (LambdaInterface) potential;
    logger.info(sb.toString());

    int nSnapshots = openers[0].countNumModels();

    if (histogramName.isEmpty()) {
      histogramName = FilenameUtils.removeExtension(filenames.get(0)) + ".his";
    }

    if (lamBins < 1) {
      File histogramFile = new File(histogramName);
      HistogramData histogramData = HistogramData.readHistogram(histogramFile);
      lamBins = histogramData.getLambdaBins();
      if (histogramData.wasHistogramRead()) {
        logger.info(format(" Autodetected %d from histogram file.", lamBins));
      }
    }

    energiesL = new ArrayList<>(lamBins);
    energiesUp = new ArrayList<>(lamBins);
    energiesDown = new ArrayList<>(lamBins);
    for (int i = 0; i < lamBins; i++) {
      energiesL.add(new ArrayList<Double>());
      energiesUp.add(new ArrayList<Double>());
      energiesDown.add(new ArrayList<Double>());
    }

    lamSep = 1.0 / (lamBins - 1);
    halfLamSep = 0.5 * lamSep;
    lamPoints = new double[lamBins];
    // TODO: Remove assumption that it's using discrete lambda bins.
    for (int i = 0; i < (lamBins - 1); i++) {
      lamPoints[i] = i * lamSep;
    }
    lamPoints[lamBins - 1] = 1.0; // Eliminate machine precision error.

    OptionalDouble optLam = openers[0].getLastReadLambda();
    // Note: OptionalDouble.isEmpty() is a JDK 11 feature, so !OptionalDouble.isPresent() preserves JDK 8 compatibility.
    if (!optLam.isPresent()) {
      throw new IllegalArgumentException(
          format(" No lambda records found in the first header of archive file %s", filenames.get(0)));
    }

    start = startFrame - 1;
    if (finalFrame < 1) {
      end = nSnapshots;
    } else {
      end = min(nSnapshots, finalFrame);
    }
    end -= startFrame; // Will always be compared to an index offset by start.

    double lambda = optLam.getAsDouble();
    int nVar = potential.getNumberOfVariables();
    x = new double[nVar];

    // --lambdaSorted and the observations array are there largely to deal with a test case that was just regular BAR with concatenated .arc files.
    observations = new int[lamBins];
    if (lambdaSorted) {
      fill(observations, -startFrame);
    } else {
      fill(observations, 0);
    }

    logger.info(" Reading snapshots.");

    addEntries(lambda, 0);

    for (int i = 1; i < end; i++) {
      for (int j = 0; j < numTopologies; j++) {
        openers[j].readNext(false, false);
      }
      lambda = openers[0].getLastReadLambda().getAsDouble();
      addEntries(lambda, i);
    }

    for (SystemFilter opener : openers) {
      opener.closeReader();
    }

    double[][] eLow = new double[lamBins][];
    double[][] eAt = new double[lamBins][];
    double[][] eHigh = new double[lamBins][];
    for (int i = 0; i < lamBins; i++) {
      eLow[i] = energiesDown.get(i).stream().mapToDouble(Double::doubleValue).toArray();
      eAt[i] = energiesL.get(i).stream().mapToDouble(Double::doubleValue).toArray();
      eHigh[i] = energiesUp.get(i).stream().mapToDouble(Double::doubleValue).toArray();
    }

    logger.info("\n Initial estimate via the iteration method.");
    BennettAcceptanceRatio bar = new BennettAcceptanceRatio(lamPoints, eLow, eAt, eHigh,
        new double[]{temp});
    Zwanzig forwards = bar.getInitialForwardsGuess();
    Zwanzig backwards = bar.getInitialBackwardsGuess();

    logger.info(format(" Free energy via BAR:           %15.9f +/- %.9f kcal/mol.",
        bar.getTotalFreeEnergyDifference(), bar.getTotalFEDifferenceUncertainty()));
    logger.info(format(" Free energy via forwards FEP:  %15.9f +/- %.9f kcal/mol.",
        forwards.getTotalFreeEnergyDifference(), forwards.getTotalFEDifferenceUncertainty()));
    logger.info(format(" Free energy via backwards FEP: %15.9f +/- %.9f kcal/mol.",
        backwards.getTotalFreeEnergyDifference(), backwards.getTotalFEDifferenceUncertainty()));
    logger.info(" Note - non-bootstrap FEP uncertainties are currently unreliable.");

    double[] barFE = bar.getFreeEnergyDifferences();
    double[] barVar = bar.getFEDifferenceUncertainties();
    double[] forwardsFE = forwards.getFreeEnergyDifferences();
    double[] forwardsVar = forwards.getFEDifferenceUncertainties();
    double[] backwardsFE = backwards.getFreeEnergyDifferences();
    double[] backwardsVar = backwards.getFEDifferenceUncertainties();

    sb = new StringBuilder(
        "\n Free Energy Profile Per Window\n Min_Lambda Counts Max_Lambda Counts         BAR_dG      BAR_Var          FEP_dG      FEP_Var     FEP_Back_dG FEP_Back_Var\n");
    for (int i = 0; i < (lamBins - 1); i++) {
      sb.append(format(" %-10.8f %6d %-10.8f %6d %15.9f %12.9f %15.9f %12.9f %15.9f %12.9f\n",
          lamPoints[i], eAt[i].length, lamPoints[i + 1], eAt[i + 1].length, barFE[i],
          barVar[i], forwardsFE[i], forwardsVar[i], backwardsFE[i], backwardsVar[i]));
    }
    logger.info(sb.toString());

    if (bootstrap == -1) {
      int totalRead = Arrays.stream(observations).min().getAsInt();
      if (totalRead >= MIN_BOOTSTRAP_TRIALS) {
        bootstrap = AUTO_BOOTSTRAP_NUMERATOR / totalRead;
        // Weird Groovy syntax because Groovy defaults to BigDecimal/BigInteger.
        bootstrap = max(MIN_BOOTSTRAP_TRIALS, min(MAX_BOOTSTRAP_TRIALS, bootstrap));
      } else {
        logger.info(format(
            " At least one lambda window had only %d snapshots read; defaulting to %d bootstrap cycles!",
            totalRead, MIN_BOOTSTRAP_TRIALS));
        bootstrap = MIN_BOOTSTRAP_TRIALS;
      }
    }

    long bootPrint = BOOTSTRAP_PRINT;
    if (!verbose) {
      bootPrint *= 10L;
    }

    // If bootstrap <= 0, skip this section.
    if (bootstrap > 0) {
      logger.info(format(" Re-estimate free energy and uncertainty from %d bootstrap trials.", bootstrap));

      EstimateBootstrapper barBS = new EstimateBootstrapper(bar);
      EstimateBootstrapper forBS = new EstimateBootstrapper(forwards);
      EstimateBootstrapper backBS = new EstimateBootstrapper(backwards);

      long time = -System.nanoTime();
      barBS.bootstrap(bootstrap, bootPrint);
      time += System.nanoTime();
      logger.info(format(" BAR bootstrapping complete in %.4f sec", time * Constants.NS2SEC));

      time = -System.nanoTime();
      forBS.bootstrap(bootstrap, bootPrint);
      time += System.nanoTime();
      logger.info(format(" Forwards FEP bootstrapping complete in %.4f sec", time * Constants.NS2SEC));

      time = -System.nanoTime();
      backBS.bootstrap(bootstrap, bootPrint);
      time += System.nanoTime();
      logger.info(format(" Reverse FEP bootstrapping complete in %.4f sec", time * Constants.NS2SEC));

      barFE = barBS.getFreeEnergyDifferences();
      barVar = barBS.getFEDifferenceStdDevs();
      forwardsFE = forBS.getFreeEnergyDifferences();
      forwardsVar = forBS.getFEDifferenceStdDevs();
      backwardsFE = backBS.getFreeEnergyDifferences();
      backwardsVar = backBS.getFEDifferenceStdDevs();

      double sumFE = barBS.getTotalFreeEnergyDifference(barFE);
      double varFE = barBS.getTotalFEDifferenceUncertainty(barVar);
      logger.info(format(" Free energy via BAR:           %15.9f +/- %.9f kcal/mol.", sumFE, varFE));

      sumFE = forBS.getTotalFreeEnergyDifference(forwardsFE);
      varFE = forBS.getTotalFEDifferenceUncertainty();
      logger.info(format(" Free energy via forwards FEP:  %15.9f +/- %.9f kcal/mol.", sumFE, varFE));

      sumFE = backBS.getTotalFreeEnergyDifference(backwardsFE);
      varFE = backBS.getTotalFEDifferenceUncertainty();
      logger.info(format(" Free energy via backwards FEP:  %15.9f +/- %.9f kcal/mol.", sumFE, varFE));

      sb = new StringBuilder(
          " Free Energy Profile\n Min_Lambda Counts Max_Lambda Counts         BAR_dG      BAR_Var          FEP_dG      FEP_Var     FEP_Back_dG FEP_Back_Var\n");
      for (int i = 0; i < (lamBins - 1); i++) {
        sb.append(format(" %-10.8f %6d %-10.8f %6d %15.9f %12.9f %15.9f %12.9f %15.9f %12.9f\n",
            lamPoints[i], eAt[i].length, lamPoints[i + 1], eAt[i + 1].length, barFE[i], barVar[i],
            forwardsFE[i], forwardsVar[i], backwardsFE[i], backwardsVar[i]));
      }
      logger.info(sb.toString());
    } else {
      logger.info(" Bootstrap resampling disabled.");
    }

    return this;
  }

  /**
   * Adds entries to the energy lists to be sent to the statistical estimators. This is where -s, --fi,
   * --st and --lambdaSorted are applied.
   *
   * @param lambda Lambda of the snapshot just read.
   * @param index  Index of this snapshot in the entire archive.
   */
  private void addEntries(double lambda, int index) {
    int bin = binForLambda(lambda);
    ++observations[bin];
    // The observation count is pre-offset by -start.
    int offsetIndex = lambdaSorted ? observations[bin] : index - start;
    assert offsetIndex <= end;

    boolean inRange = offsetIndex >= 0 && offsetIndex <= end;
    boolean onStride = (offsetIndex % stride == 0);
    if (inRange && onStride) {
      x = potential.getCoordinates(x);
      lastEntries[0] = addLambdaDown(lambda, bin);
      lastEntries[1] = addAtLambda(lambda, bin);
      lastEntries[2] = addLambdaUp(lambda, bin);

      String low = Double.isNaN(lastEntries[0]) ? nanFormat : format(energyFormat, lastEntries[0]);
      String high = Double.isNaN(lastEntries[2]) ? nanFormat : format(energyFormat, lastEntries[2]);

      logger.log(standardLogging, format(" Energies for snapshot %5d at lambda %.4f: " +
          "%s, %s, %s", (index + 1), lambda, low, format(energyFormat, lastEntries[1]), high));
    } else {
      logger.log(standardLogging, " Skipping frame " + index);
    }
  }

  /**
   * Adds an entry to the energiesL list.
   *
   * @param lambda Lambda of the last read snapshot.
   * @param bin    Lambda bin of this snapshot.
   * @return Energy at lambda = lambda.
   */
  private double addAtLambda(double lambda, int bin) {
    assert lambda >= 0.0 && lambda <= 1.0;
    linter.setLambda(lambda);
    double e = potential.energy(x, false);
    // Following line is only for debugging purposes!
    //e += DEBUG_OFFSET * lambda;
    energiesL.get(bin).add(e);
    return e;
  }

  /**
   * Adds an entry to the energiesUp list.
   *
   * @param lambda Lambda of the last read snapshot.
   * @param bin    Lambda bin of this snapshot.
   * @return Energy at lambda = lambda+dL.
   */
  private double addLambdaUp(double lambda, int bin) {
    double modLambda = lambda + lamSep;
    // DISCRETE ONLY: assert lambda == 1.0d || modLambda < (1.0 + 1.0E-6);
    modLambda = min(1.0d, modLambda);
    if (bin == (lamBins - 1)) {
      energiesUp.get(bin).add(Double.NaN);
      return Double.NaN;
    } else {
      linter.setLambda(modLambda);
      double e = potential.energy(x, false);
      // Following line is only for debugging purposes!
      //e += DEBUG_OFFSET * modLambda;
      energiesUp.get(bin).add(e);
      linter.setLambda(lambda);
      return e;
    }
  }

  /**
   * Adds an entry to the energiesDown list.
   *
   * @param lambda Lambda of the last read snapshot.
   * @param bin    Lambda bin of this snapshot.
   * @return Energy at lambda = lambda-dL.
   */
  private double addLambdaDown(double lambda, int bin) {
    double modLambda = lambda - lamSep;
    // DISCRETE ONLY: assert lambda == 0.0d || modLambda > -1.0E-6;
    modLambda = max(0.0d, modLambda);

    if (bin == 0) {
      energiesDown.get(0).add(Double.NaN);
      return Double.NaN;
    } else {
      linter.setLambda(modLambda);
      double e = potential.energy(x, false);
      // Following line is only for debugging purposes!
      //e += DEBUG_OFFSET * modLambda;
      energiesDown.get(bin).add(e);
      linter.setLambda(lambda);
      return e;
    }
  }

  /**
   * <p>binForLambda.</p>
   *
   * @param lambda a double.
   * @return a int.
   */
  private int binForLambda(double lambda) {
    // TODO: Robust to existence of half-bins.
    return (int) round(lambda / lamSep);
  }
}
