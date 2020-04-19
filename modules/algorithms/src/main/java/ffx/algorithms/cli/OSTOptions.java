// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.algorithms.cli;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.thermodynamics.HistogramSettings;
import ffx.algorithms.thermodynamics.MonteCarloOST;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering.OptimizationParameters;
import ffx.crystal.CrystalPotential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.cli.WriteoutOptions;
import ffx.utilities.Constants;
import java.io.File;
import java.io.IOException;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.configuration2.Configuration;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that utilize variants of the Orthogonal Space
 * Tempering (OST) algorithm. Metadynamics will be treated as a special case of OST where there is
 * no dU/dL axis.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class OSTOptions {

  private static final Logger logger = Logger.getLogger(OSTOptions.class.getName());

  /**
   * -c or --count Sets the number of time steps between OST counts.
   */
  @Option(
      names = {"-C", "--count"},
      paramLabel = "10",
      defaultValue = "10",
      description = "Time steps between MD Orthogonal Space counts.")
  public int countFreq;

  /**
   * --bM or --biasMag sets the initial Gaussian bias magnitude in kcal/mol.
   */
  @Option(
      names = {"--bM", "--biasMag"},
      paramLabel = "0.05",
      defaultValue = "0.05",
      split = ",",
      description =
          "Orthogonal Space Gaussian bias magnitude (kcal/mol); RepEx OST uses a comma-separated list.")
  private double[] biasMag;

  /**
   * --iW or --independentWalkers enforces that each walker maintains their own histogram.
   */
  @Option(
      names = {"--iW", "--independentWalkers"},
      defaultValue = "false",
      description = "Enforces that each walker maintains their own histogram. ")
  private boolean independentWalkers;

  /**
   * --tp or --temperingParam sets the Dama et al tempering rate parameter, in multiples of kBT.
   */
  @Option(
      names = {"--tp", "--temperingParam"},
      paramLabel = "4.0",
      defaultValue = "4.0",
      split = ",",
      description =
          "Tempering rate parameter in multiples of kBT; repex OST uses a comma-separated list.")
  private double[] temperParam;

  /**
   * --tth or --temperingThreshold sets the tempering threshold/offset in kcal/mol.
   */
  @Option(
      names = {"--tth", "--temperingThreshold"},
      paramLabel = "20*bias",
      defaultValue = "-1",
      split = ",",
      description = "Tempering threshold in kcal/mol; repex OST uses a comma-separated list.")
  private double[] temperThreshold;

  /**
   * --mc or --monteCarlo sets the Monte Carlo scheme for Orthogonal Space Tempering.
   */
  @Option(
      names = {"--mc", "--monteCarlo"},
      defaultValue = "false",
      description = "Specify use of Monte Carlo OST")
  private boolean mc;

  /**
   * --mcHW or --monteCarloHardWall sets the Monte Carlo scheme to use a hard wall that rejects any
   * sample (Lambda, dU/dL) located in an empty histogram bin.
   */
  @Option(
      names = {"--mcHW", "--monteCarloHardWall"},
      defaultValue = "false",
      description = "Monte Carlo OST hard wall constraint.")
  private boolean mcHW;

  /**
   * --mcmD or --mcTraj Sets the number of steps to take for each MD trajectory for MC-OST.
   */
  @Option(
      names = {"--mcMD", "--mcTraj"},
      paramLabel = "100",
      defaultValue = "100",
      description = "Number of dynamics steps to take for each MD trajectory for Monte Carlo OST")
  private int mcMD;

  /**
   * --mcL or --mcLambdaStd Sets the standard deviation for lambda moves.
   */
  @Option(
      names = {"--mcL", "--mcLambdaStd"},
      paramLabel = "0.01",
      defaultValue = "0.01",
      description = "Standard deviation for lambda move.")
  private double mcL;

  /**
   * --ts or --twoStep MC Orthogonal Space sampling using separate lambda and MD moves.
   */
  @Option(
      names = {"--ts", "--twoStep"},
      defaultValue = "false",
      description = "MC Orthogonal Space sampling using separate lambda and MD moves.")
  private boolean ts;

  /**
   * --lw or --lambdaWriteOut Only write out snapshots if lambda is greater than the value
   * specified.
   */
  @Option(
      names = {"--lw", "--lambdaWriteOut"},
      paramLabel = "0.0",
      defaultValue = "0.0",
      description = "Only write out snapshots if lambda is greater than the value specified.")
  private double lambdaWriteOut;

  /**
   * Static method for constructing an orthogonal space tempering object with numerous default
   * settings. Largely indended for use with the Histogram script and other scripts which don't need
   * to actively perform OST, just read its histogram.
   *
   * @param potential a {@link ffx.crystal.CrystalPotential} to build the OST around.
   * @param lambdaRestart a {@link java.io.File} lambda restart file.
   * @param histogramRestart a {@link java.io.File} histogram restart file.
   * @param firstAssembly the first {@link ffx.potential.MolecularAssembly} in the OST system.
   * @param addedProperties a {@link org.apache.commons.configuration2.Configuration} with
   *     additional properties.
   * @param aListener any {@link ffx.algorithms.AlgorithmListener} that OST should update.
   * @return the newly built {@link OrthogonalSpaceTempering} object.
   * @throws IOException Can be thrown by errors reading restart files.
   */
  public static OrthogonalSpaceTempering constructOST(
      CrystalPotential potential,
      File lambdaRestart,
      File histogramRestart,
      MolecularAssembly firstAssembly,
      Configuration addedProperties,
      AlgorithmListener aListener)
      throws IOException {
    LambdaInterface linter = (LambdaInterface) potential;
    CompositeConfiguration allProperties =
        new CompositeConfiguration(firstAssembly.getProperties());
    if (addedProperties != null) {
      allProperties.addConfiguration(addedProperties);
    }

    String lamRestartName =
        lambdaRestart == null
            ? histogramRestart.toString().replaceFirst("\\.his$", ".lam")
            : lambdaRestart.toString();
    HistogramSettings hOps = new HistogramSettings(histogramRestart, lamRestartName, allProperties);

    // These fields are needed for the OST constructor, but otherwise are not used.
    boolean asynchronous = false;
    double timeStep = 1.0;
    double printInterval = 1.0;
    double saveInterval = 100.0;
    double temperature = 298.15;

    return new OrthogonalSpaceTempering(
        linter,
        potential,
        lambdaRestart,
        hOps,
        allProperties,
        temperature,
        timeStep,
        printInterval,
        saveInterval,
        asynchronous,
        false,
        aListener,
        0.0);
  }

  /**
   * Applies relevant options to a OST, and returns either the OST object or something that wraps
   * the OST (such as a Barostat).
   *
   * @param orthogonalSpaceTempering Orthogonal Space Tempering.
   * @param firstAssembly Primary assembly in OST.
   * @param dynamics MD options.
   * @param barostat NPT options.
   * @return a {@link ffx.crystal.CrystalPotential} object.
   */
  public CrystalPotential applyAllOSTOptions(
      OrthogonalSpaceTempering orthogonalSpaceTempering,
      MolecularAssembly firstAssembly,
      DynamicsOptions dynamics,
      BarostatOptions barostat) {
    if (dynamics.getOptimize()) {
      OptimizationParameters opt = orthogonalSpaceTempering.getOptimizationParameters();
      opt.setOptimization(true, firstAssembly);
    }
    return barostat.checkNPT(firstAssembly, orthogonalSpaceTempering);
  }

  /**
   * Assembles a MolecularDynamics wrapped around a Potential.
   *
   * @param topologies MolecularAssembly[]
   * @param potential Potential to run on
   * @param dynamics DynamicsOptions
   * @param aListener AlgorithmListener
   * @return MolecularDynamics
   */
  public MolecularDynamics assembleMolecularDynamics(
      MolecularAssembly[] topologies,
      CrystalPotential potential,
      DynamicsOptions dynamics,
      AlgorithmListener aListener) {
    // Create the MolecularDynamics instance.
    MolecularAssembly firstTop = topologies[0];
    CompositeConfiguration props = firstTop.getProperties();

    dynamics.init();

    MolecularDynamics molDyn =
        MolecularDynamics.dynamicsFactory(
            firstTop,
            potential,
            props,
            aListener,
            dynamics.thermostat,
            dynamics.integrator,
            MolecularDynamics.DynamicsEngine.FFX);
    for (int i = 1; i < topologies.length; i++) {
      molDyn.addAssembly(topologies[i], topologies[i].getProperties());
    }
    molDyn.setRestartFrequency(dynamics.getCheckpoint());

    return molDyn;
  }

  /**
   * G Runs MC-OST.
   *
   * @param monteCarloOST         MC-OST to run.
   * @param dynamicsOptions       Dynamics options.
   * @param thermodynamicsOptions Thermodynamics options.
   */
  public void beginMCOST(
      MonteCarloOST monteCarloOST,
      DynamicsOptions dynamicsOptions,
      ThermodynamicsOptions thermodynamicsOptions) {
    long nEquil = thermodynamicsOptions.getEquilSteps();

    if (nEquil > 0) {
      logger.info("\n Beginning MC-OST equilibration.");
      monteCarloOST.setMDMoveParameters(nEquil);
      if (ts) {
        monteCarloOST.sampleTwoStep();
      } else {
        monteCarloOST.sampleOneStep();
      }
      monteCarloOST.setEquilibration(false);
      logger.info("\n Finished MC-OST equilibration.");
    }

    logger.info("\n Beginning MC-OST sampling.");
    monteCarloOST.setLambdaStdDev(mcL);
    monteCarloOST.setMDMoveParameters(dynamicsOptions.steps);

    if (ts) {
      monteCarloOST.sampleTwoStep();
    } else {
      monteCarloOST.sampleOneStep();
    }
  }

  /**
   * Begins MD-OST sampling from an assembled OST.
   *
   * @param orthogonalSpaceTempering The OST object.
   * @param molecularAssemblies      All MolecularAssemblys.
   * @param crystalPotential         The top-layer CrystalPotential.
   * @param dynamicsOptions          Dynamics options.
   * @param writeOut                 a {@link WriteoutOptions} object.
   * @param thermo                   Thermodynamics options.
   * @param dynFile                  The .dyn dynamics restart file.
   * @param algorithmListener        AlgorithmListener
   */
  public void beginMDOST(
      OrthogonalSpaceTempering orthogonalSpaceTempering,
      MolecularAssembly[] molecularAssemblies,
      CrystalPotential crystalPotential,
      DynamicsOptions dynamicsOptions,
      WriteoutOptions writeOut,
      ThermodynamicsOptions thermo,
      File dynFile,
      AlgorithmListener algorithmListener) {

    dynamicsOptions.init();

    MolecularDynamics molDyn =
        assembleMolecularDynamics(
            molecularAssemblies, crystalPotential, dynamicsOptions, algorithmListener);

    boolean initVelocities = true;
    long nSteps = dynamicsOptions.steps;
    // Start sampling.
    long nEquil = thermo.getEquilSteps();
    if (nEquil > 0) {
      logger.info("\n Beginning equilibration");
      orthogonalSpaceTempering.setPropagateLambda(false);
      runDynamics(molDyn, nEquil, dynamicsOptions, writeOut, true, dynFile);
      logger.info(" Beginning OST sampling");
      orthogonalSpaceTempering.setPropagateLambda(true);
    } else {
      logger.info(" Beginning OST sampling without equilibration");
      if (!thermo.getResetNumSteps()) {
        long nEnergyCount = orthogonalSpaceTempering.getEnergyCount();
        if (nEnergyCount > 0) {
          nSteps -= nEnergyCount;
          logger.info(
              String.format(
                  " Lambda file: %12d steps picked up, now sampling %12d steps",
                  nEnergyCount, nSteps));
          initVelocities = false;
        }
      }
    }
    if (nSteps > 0) {
      runDynamics(molDyn, nSteps, dynamicsOptions, writeOut, initVelocities, dynFile);
    } else {
      logger.info(" No steps remaining for this process!");
    }
  }

  /**
   * constructOST.
   *
   * @param crystalPotential a {@link ffx.crystal.CrystalPotential} to build the OST around.
   * @param lambdaRestartFile a {@link java.io.File} lambda restart file.
   * @param histogramRestartFile a {@link java.io.File} histogram restart file.
   * @param firstAssembly the first {@link ffx.potential.MolecularAssembly} in the OST system.
   * @param addedProperties a {@link org.apache.commons.configuration2.Configuration} with
   *     additional properties.
   * @param dynamicsOptions a {@link ffx.algorithms.cli.DynamicsOptions} with MD-related settings.
   * @param thermodynamicsOptions a {@link ffx.algorithms.cli.ThermodynamicsOptions} with
   *     thermodynamics-related settings.
   * @param lambdaParticleOptions a {@link ffx.algorithms.cli.LambdaParticleOptions} with lambda
   *     particle-related settings.
   * @param algorithmListener any {@link ffx.algorithms.AlgorithmListener} that OST should update.
   * @param async If OST should use asynchronous communications.
   * @return the newly built {@link OrthogonalSpaceTempering} object.
   * @throws IOException Can be thrown by errors reading restart files.
   */
  public OrthogonalSpaceTempering constructOST(
      CrystalPotential crystalPotential,
      File lambdaRestartFile,
      File histogramRestartFile,
      MolecularAssembly firstAssembly,
      Configuration addedProperties,
      DynamicsOptions dynamicsOptions,
      ThermodynamicsOptions thermodynamicsOptions,
      LambdaParticleOptions lambdaParticleOptions,
      AlgorithmListener algorithmListener,
      boolean async)
      throws IOException {

    LambdaInterface lambdaInterface = (LambdaInterface) crystalPotential;
    CompositeConfiguration compositeConfiguration =
        new CompositeConfiguration(firstAssembly.getProperties());
    if (addedProperties != null) {
      compositeConfiguration.addConfiguration(addedProperties);
    }
    double temp = dynamicsOptions.getTemp();
    double dt = dynamicsOptions.getDt();
    double report = dynamicsOptions.getReport();
    double checkpoint = dynamicsOptions.getCheckpoint();
    boolean resetNSteps = thermodynamicsOptions.getResetNumSteps();

    String lamRestartName =
        lambdaRestartFile == null
            ? histogramRestartFile.toString().replaceFirst("\\.his$", ".lam")
            : lambdaRestartFile.toString();
    HistogramSettings histogramSettings =
        generateHistogramSettings(
            histogramRestartFile,
            lamRestartName,
            compositeConfiguration,
            0,
            dynamicsOptions,
            lambdaParticleOptions,
            independentWalkers,
            async);

    OrthogonalSpaceTempering orthogonalSpaceTempering =
        new OrthogonalSpaceTempering(
            lambdaInterface,
            crystalPotential,
            lambdaRestartFile,
            histogramSettings,
            compositeConfiguration,
            temp,
            dt,
            report,
            checkpoint,
            async,
            resetNSteps,
            algorithmListener,
            lambdaWriteOut);
    orthogonalSpaceTempering.setHardWallConstraint(mcHW);

    // Do NOT run applyOSTOptions here, because that can mutate the OST to a Barostat.
    return orthogonalSpaceTempering;
  }

  /**
   * Constructs histogram settings based on stored values.
   *
   * @param histogramRestartFile Histogram restart (.his) file.
   * @param lambdaFileName Name of the lambda file.
   * @param allProperties All properties requested.
   * @param index Index of the histogram (repex/independent walkers).
   * @param dynamicsOptions Dynamics options.
   * @param lambdaParticleOptions Lambda particle options.
   * @param writeIndependent Whether walkers should write their own histogram files.
   * @param async Whether to use asynchronous communication.
   * @return Settings to apply to a new Histogram object.
   * @throws IOException From reading the histogram file.
   */
  public HistogramSettings generateHistogramSettings(
      File histogramRestartFile,
      String lambdaFileName,
      CompositeConfiguration allProperties,
      int index,
      DynamicsOptions dynamicsOptions,
      LambdaParticleOptions lambdaParticleOptions,
      boolean writeIndependent,
      boolean async)
      throws IOException {
    return generateHistogramSettings(
        histogramRestartFile,
        lambdaFileName,
        allProperties,
        index,
        dynamicsOptions,
        lambdaParticleOptions,
        writeIndependent,
        async,
        false);
  }

  /**
   * Constructs histogram settings based on stored values.
   *
   * @param histogramRestartFile Histogram restart (.his) file.
   * @param lambdaFileName Name of the lambda file.
   * @param allProperties All properties requested.
   * @param index Index of the histogram (repex/independent walkers).
   * @param dynamics Dynamics options.
   * @param lambdaParticleOptions Lambda particle options.
   * @param writeIndependent Whether walkers should write their own histogram files.
   * @param async Whether to use asynchronous communication.
   * @param overrideHistogram Whether to override settings stored in the histogram restart file (if
   *     it exists).
   * @return Settings to apply to a new Histogram object.
   * @throws IOException From reading the histogram file.
   */
  public HistogramSettings generateHistogramSettings(
      File histogramRestartFile,
      String lambdaFileName,
      CompositeConfiguration allProperties,
      int index,
      DynamicsOptions dynamics,
      LambdaParticleOptions lambdaParticleOptions,
      boolean writeIndependent,
      boolean async,
      boolean overrideHistogram)
      throws IOException {
    HistogramSettings histogramSettings =
        new HistogramSettings(histogramRestartFile, lambdaFileName, allProperties);
    histogramSettings.temperingFactor = getTemperingParameter(index);
    if (thresholdsSet()) {
      histogramSettings.setTemperOffset(getTemperingThreshold(index));
    }
    histogramSettings.dt = dynamics.getDt() * Constants.FSEC_TO_PSEC;
    histogramSettings.setWriteIndependent(writeIndependent);
    histogramSettings.setIndependentWalkers(independentWalkers);
    histogramSettings.asynchronous = async;
    if (overrideHistogram || !histogramRestartFile.exists()) {
      histogramSettings.setBiasMag(getBiasMag(index));
      // TODO: Let temperature be an array thing too.
      histogramSettings.temperature = dynamics.getTemp();
      histogramSettings.thetaFriction = lambdaParticleOptions.getLambdaFriction();
      histogramSettings.thetaMass = lambdaParticleOptions.getLambdaMass();
      histogramSettings.countInterval = countFreq;
    }

    return histogramSettings;
  }

  /**
   * Checks if independent walkers has been specified.
   *
   * @return Walker independence.
   */
  public boolean getIndependentWalkers() {
    return independentWalkers;
  }

  /**
   * Checks if use of the Monte Carlo algorithm has been specified.
   *
   * @return Monte Carlo OST (as opposed to molecular dynamics OST).
   */
  public boolean isMc() {
    return mc;
  }

  /**
   * Returns true if the 2-step option is enabled (not guaranteed to also mean that MC is enabled!).
   *
   * @return If --ts is enabled.
   */
  public boolean isTwoStep() {
    return ts;
  }

  /**
   * setupMCOST.
   *
   * @param orthogonalSpaceTempering a {@link OrthogonalSpaceTempering} object.
   * @param molecularAssemblies      an array of {@link ffx.potential.MolecularAssembly} objects.
   * @param dynamicsOptions          a {@link ffx.algorithms.cli.DynamicsOptions} object.
   * @param thermodynamicsOptions    a {@link ffx.algorithms.cli.ThermodynamicsOptions} object.
   * @param verbose                  Whether to print out additional information about MC-OST.
   * @param algorithmListener        An AlgorithmListener
   * @return An assembled MonteCarloOST ready to run.
   */
  public MonteCarloOST setupMCOST(
      OrthogonalSpaceTempering orthogonalSpaceTempering,
      MolecularAssembly[] molecularAssemblies,
      DynamicsOptions dynamicsOptions,
      ThermodynamicsOptions thermodynamicsOptions,
      boolean verbose,
      AlgorithmListener algorithmListener) {
    dynamicsOptions.init();

    MonteCarloOST monteCarloOST =
        new MonteCarloOST(
            orthogonalSpaceTempering.getPotentialEnergy(),
            orthogonalSpaceTempering,
            molecularAssemblies[0],
            molecularAssemblies[0].getProperties(),
            algorithmListener,
            dynamicsOptions,
            verbose,
            mcMD);

    MolecularDynamics md = monteCarloOST.getMD();
    for (int i = 1; i < molecularAssemblies.length; i++) {
      md.addAssembly(molecularAssemblies[i], molecularAssemblies[i].getProperties());
    }

    long nEquil = thermodynamicsOptions.getEquilSteps();
    if (nEquil > 0) {
      monteCarloOST.setEquilibration(true);
    }
    return monteCarloOST;
  }

  private boolean thresholdsSet() {
    return temperThreshold.length != 1 || temperThreshold[0] >= 0;
  }

  private void runDynamics(
      MolecularDynamics molecularDynamics,
      long numSteps,
      DynamicsOptions dynamicsOptions,
      WriteoutOptions writeoutOptions,
      boolean initVelocities,
      File dyn) {
    molecularDynamics.dynamic(
        numSteps,
        dynamicsOptions.dt,
        dynamicsOptions.report,
        dynamicsOptions.write,
        dynamicsOptions.temp,
        initVelocities,
        writeoutOptions.getFileType(),
        dynamicsOptions.getCheckpoint(),
        dyn);
  }

  /**
   * Returns the initial bias magnitude associated with a walker.
   *
   * @param i Index of a walker
   * @return Its intended initial bias magnitude in kcal/mol.
   */
  private double getBiasMag(int i) {
    return biasMag.length > 1 ? biasMag[i] : biasMag[0];
  }

  /**
   * Returns the tempering threshold associated with a walker.
   *
   * @param i Index of a walker
   * @return Its intended tempering threshold in kcal/mol.
   */
  private double getTemperingThreshold(int i) {
    return temperThreshold.length > 1 ? temperThreshold[i] : temperThreshold[0];
  }

  /**
   * Returns the tempering parameter associated with a walker.
   *
   * @param i Index of a walker
   * @return Its intended tempering parameter in kBT.
   */
  private double getTemperingParameter(int i) {
    return temperParam.length > 1 ? temperParam[i] : temperParam[0];
  }
}
