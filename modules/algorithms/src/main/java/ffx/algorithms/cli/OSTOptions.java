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
import picocli.CommandLine.ArgGroup;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that utilize variants of the Orthogonal Space
 * Tempering (OST) algorithm. Metadynamics will be treated as a special case of OST where there is no
 * dU/dL axis.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class OSTOptions {

  private static final Logger logger = Logger.getLogger(OSTOptions.class.getName());

  /**
   * The ArgGroup keeps the OSTOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Orthogonal Space Tempering Options%n", validate = false)
  public OSTOptionGroup group = new OSTOptionGroup();

  /**
   * The ArgGroup keeps the OSTOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Monte Carlo Orthogonal Space Tempering Options%n", validate = false)
  public MCOSTOptionGroup mcGroup = new MCOSTOptionGroup();

  /**
   * Static method for constructing an orthogonal space tempering object with numerous default
   * settings. Largely indended for use with the Histogram script and other scripts which don't need
   * to actively perform OST, just read its histogram.
   *
   * @param crystalPotential a {@link ffx.crystal.CrystalPotential} to build the OST around.
   * @param lambdaRestart a {@link java.io.File} lambda restart file.
   * @param histogramRestart a {@link java.io.File} histogram restart file.
   * @param firstAssembly the first {@link ffx.potential.MolecularAssembly} in the OST system.
   * @param configuration a {@link org.apache.commons.configuration2.Configuration} with
   *     additional properties.
   * @param algorithmListener any {@link ffx.algorithms.AlgorithmListener} that OST should
   *     update.
   * @return the newly built {@link OrthogonalSpaceTempering} object.
   * @throws IOException Can be thrown by errors reading restart files.
   */
  public static OrthogonalSpaceTempering constructOST(
      CrystalPotential crystalPotential,
      File lambdaRestart,
      File histogramRestart,
      MolecularAssembly firstAssembly,
      Configuration configuration,
      AlgorithmListener algorithmListener)
      throws IOException {
    LambdaInterface linter = (LambdaInterface) crystalPotential;
    CompositeConfiguration allProperties =
        new CompositeConfiguration(firstAssembly.getProperties());
    if (configuration != null) {
      allProperties.addConfiguration(configuration);
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
        crystalPotential,
        lambdaRestart,
        hOps,
        allProperties,
        temperature,
        timeStep,
        printInterval,
        saveInterval,
        asynchronous,
        false,
        algorithmListener,
        0.0);
  }

  /**
   * Applies relevant options to a OST, and returns either the OST object or something that wraps the
   * OST (such as a Barostat).
   *
   * @param orthogonalSpaceTempering Orthogonal Space Tempering.
   * @param firstAssembly Primary assembly in OST.
   * @param dynamicsOptions MD options.
   * @param barostatOptions NPT options.
   * @return a {@link ffx.crystal.CrystalPotential} object.
   */
  public CrystalPotential applyAllOSTOptions(
      OrthogonalSpaceTempering orthogonalSpaceTempering,
      MolecularAssembly firstAssembly,
      DynamicsOptions dynamicsOptions,
      BarostatOptions barostatOptions) {
    if (dynamicsOptions.getOptimize()) {
      OptimizationParameters opt = orthogonalSpaceTempering.getOptimizationParameters();
      opt.setOptimization(true, firstAssembly);
    }
    return barostatOptions.checkNPT(firstAssembly, orthogonalSpaceTempering);
  }

  /**
   * G Runs MC-OST.
   *
   * @param monteCarloOST MC-OST to run.
   * @param dynamicsOptions Dynamics options.
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
      if (mcGroup.twoStep) {
        monteCarloOST.sampleTwoStep();
      } else {
        monteCarloOST.sampleOneStep();
      }
      monteCarloOST.setEquilibration(false);
      logger.info("\n Finished MC-OST equilibration.");
    }

    logger.info("\n Beginning MC-OST sampling.");
    monteCarloOST.setLambdaStdDev(mcGroup.mcLambdaStdDev);
    monteCarloOST.setMDMoveParameters(dynamicsOptions.getSteps());

    if (mcGroup.twoStep) {
      monteCarloOST.sampleTwoStep();
    } else {
      monteCarloOST.sampleOneStep();
    }
  }

  /**
   * Assembles a MolecularDynamics wrapped around a Potential.
   *
   * @param molecularAssemblies MolecularAssembly[]
   * @param crystalPotential Potential to run on
   * @param dynamicsOptions DynamicsOptions
   * @param algorithmListener AlgorithmListener
   * @return MolecularDynamics
   */
  public MolecularDynamics assembleMolecularDynamics(
      MolecularAssembly[] molecularAssemblies,
      CrystalPotential crystalPotential,
      DynamicsOptions dynamicsOptions,
      AlgorithmListener algorithmListener) {
    // Create the MolecularDynamics instance.
    MolecularAssembly firstTop = molecularAssemblies[0];
    CompositeConfiguration props = firstTop.getProperties();

    dynamicsOptions.init();

    MolecularDynamics molDyn =
        MolecularDynamics.dynamicsFactory(
            firstTop,
            crystalPotential,
            props,
            algorithmListener,
            dynamicsOptions.thermostat,
            dynamicsOptions.integrator,
            MolecularDynamics.DynamicsEngine.FFX);
    for (int i = 1; i < molecularAssemblies.length; i++) {
      molDyn.addAssembly(molecularAssemblies[i], molecularAssemblies[i].getProperties());
    }
    molDyn.setRestartFrequency(dynamicsOptions.getCheckpoint());

    return molDyn;
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
   * @param dynamicsOptions a {@link ffx.algorithms.cli.DynamicsOptions} with MD-related
   *     settings.
   * @param thermodynamicsOptions a {@link ffx.algorithms.cli.ThermodynamicsOptions} with
   *     thermodynamics-related settings.
   * @param lambdaParticleOptions a {@link ffx.algorithms.cli.LambdaParticleOptions} with lambda
   *     particle-related settings.
   * @param algorithmListener any {@link ffx.algorithms.AlgorithmListener} that OST should
   *     update.
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
    double temp = dynamicsOptions.getTemperature();
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
            group.independentWalkers,
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
            group.lambdaWriteOut);
    orthogonalSpaceTempering.setHardWallConstraint(mcGroup.mcHardWall);

    // Do NOT run applyOSTOptions here, because that can mutate the OST to a Barostat.
    return orthogonalSpaceTempering;
  }

  /**
   * Begins MD-OST sampling from an assembled OST.
   *
   * @param orthogonalSpaceTempering The OST object.
   * @param molecularAssemblies All MolecularAssemblys.
   * @param crystalPotential The top-layer CrystalPotential.
   * @param dynamicsOptions Dynamics options.
   * @param writeoutOptions a {@link WriteoutOptions} object.
   * @param thermodynamicsOptions Thermodynamics options.
   * @param dynFile The .dyn dynamics restart file.
   * @param algorithmListener AlgorithmListener
   */
  public void beginMDOST(
      OrthogonalSpaceTempering orthogonalSpaceTempering,
      MolecularAssembly[] molecularAssemblies,
      CrystalPotential crystalPotential,
      DynamicsOptions dynamicsOptions,
      WriteoutOptions writeoutOptions,
      ThermodynamicsOptions thermodynamicsOptions,
      File dynFile,
      AlgorithmListener algorithmListener) {

    dynamicsOptions.init();

    MolecularDynamics molDyn =
        assembleMolecularDynamics(
            molecularAssemblies, crystalPotential, dynamicsOptions, algorithmListener);

    boolean initVelocities = true;
    long nSteps = dynamicsOptions.getSteps();
    // Start sampling.
    long nEquil = thermodynamicsOptions.getEquilSteps();
    if (nEquil > 0) {
      logger.info("\n Beginning equilibration");
      orthogonalSpaceTempering.setPropagateLambda(false);
      runDynamics(molDyn, nEquil, dynamicsOptions, writeoutOptions, true, dynFile);
      logger.info(" Beginning OST sampling");
      orthogonalSpaceTempering.setPropagateLambda(true);
    } else {
      logger.info(" Beginning OST sampling without equilibration");
      if (!thermodynamicsOptions.getResetNumSteps()) {
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
      runDynamics(molDyn, nSteps, dynamicsOptions, writeoutOptions, initVelocities, dynFile);
    } else {
      logger.info(" No steps remaining for this process!");
    }
  }

  /**
   * Constructs histogram settings based on stored values.
   *
   * @param histogramRestartFile Histogram restart (.his) file.
   * @param lambdaFileName Name of the lambda file.
   * @param compositeConfiguration All properties requested.
   * @param index Index of the histogram (RepEx/independent walkers).
   * @param dynamicsOptions Dynamics options.
   * @param lambdaParticleOptions Lambda particle options.
   * @param writeIndependent Whether walkers should write their own histogram files.
   * @param async Whether to use asynchronous communication.
   * @param overrideHistogram Whether to override settings stored in the histogram restart file
   *     (if it exists).
   * @return Settings to apply to a new Histogram object.
   * @throws IOException From reading the histogram file.
   */
  public HistogramSettings generateHistogramSettings(
      File histogramRestartFile,
      String lambdaFileName,
      CompositeConfiguration compositeConfiguration,
      int index,
      DynamicsOptions dynamicsOptions,
      LambdaParticleOptions lambdaParticleOptions,
      boolean writeIndependent,
      boolean async,
      boolean overrideHistogram)
      throws IOException {
    HistogramSettings histogramSettings =
        new HistogramSettings(histogramRestartFile, lambdaFileName, compositeConfiguration);
    histogramSettings.temperingFactor = getTemperingParameter(index);
    if (thresholdsSet()) {
      histogramSettings.setTemperOffset(getTemperingThreshold(index));
    }
    histogramSettings.dt = dynamicsOptions.getDt() * Constants.FSEC_TO_PSEC;
    histogramSettings.setWriteIndependent(writeIndependent);
    histogramSettings.setIndependentWalkers(group.independentWalkers);
    histogramSettings.asynchronous = async;
    if (overrideHistogram || !histogramRestartFile.exists()) {
      histogramSettings.setBiasMag(getBiasMag(index));
      // TODO: Let temperature be an array thing too.
      histogramSettings.temperature = dynamicsOptions.getTemperature();
      histogramSettings.thetaFriction = lambdaParticleOptions.getLambdaFriction();
      histogramSettings.thetaMass = lambdaParticleOptions.getLambdaMass();
      histogramSettings.countInterval = group.countInterval;
    }

    return histogramSettings;
  }

  /**
   * Constructs histogram settings based on stored values.
   *
   * @param histogramRestartFile Histogram restart (.his) file.
   * @param lambdaFileName Name of the lambda file.
   * @param compositeConfiguration All properties requested.
   * @param index Index of the histogram (RepEx/independent walkers).
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
      CompositeConfiguration compositeConfiguration,
      int index,
      DynamicsOptions dynamicsOptions,
      LambdaParticleOptions lambdaParticleOptions,
      boolean writeIndependent,
      boolean async)
      throws IOException {
    return generateHistogramSettings(
        histogramRestartFile,
        lambdaFileName,
        compositeConfiguration,
        index,
        dynamicsOptions,
        lambdaParticleOptions,
        writeIndependent,
        async,
        false);
  }

  /**
   * Checks if independent walkers has been specified.
   *
   * @return Walker independence.
   */
  public boolean getIndependentWalkers() {
    return group.independentWalkers;
  }

  /**
   * Checks if use of the Monte Carlo algorithm has been specified.
   *
   * @return Monte Carlo OST (as opposed to molecular dynamics OST).
   */
  public boolean isMonteCarlo() {
    return mcGroup.monteCarlo;
  }

  public void setMonteCarlo(boolean monteCarlo) {
    mcGroup.monteCarlo = monteCarlo;
  }

  /**
   * Returns true if the 2-step option is enabled (not guaranteed to also mean that MC is enabled!).
   *
   * @return If --ts is enabled.
   */
  public boolean isTwoStep() {
    return mcGroup.twoStep;
  }

  public void setTwoStep(boolean twoStep) {
    mcGroup.twoStep = twoStep;
  }

  /**
   * setupMCOST.
   *
   * @param orthogonalSpaceTempering a {@link OrthogonalSpaceTempering} object.
   * @param molecularAssemblies an array of {@link ffx.potential.MolecularAssembly} objects.
   * @param dynamicsOptions a {@link ffx.algorithms.cli.DynamicsOptions} object.
   * @param thermodynamicsOptions a {@link ffx.algorithms.cli.ThermodynamicsOptions} object.
   * @param verbose Whether to print out additional information about MC-OST.
   * @param algorithmListener An AlgorithmListener
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
            mcGroup.mcMDSteps);

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
    return group.temperingThreshold.length != 1 || group.temperingThreshold[0] >= 0;
  }

  /**
   * Returns the initial bias magnitude associated with a walker.
   *
   * @param i Index of a walker
   * @return Its intended initial bias magnitude in kcal/mol.
   */
  private double getBiasMag(int i) {
    return group.biasMag.length > 1 ? group.biasMag[i] : group.biasMag[0];
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
        dynamicsOptions.getDt(),
        dynamicsOptions.getReport(),
        dynamicsOptions.getWrite(),
        dynamicsOptions.getTemperature(),
        initVelocities,
        writeoutOptions.getFileType(),
        dynamicsOptions.getCheckpoint(),
        dyn);
  }

  /**
   * Returns the tempering threshold associated with a walker.
   *
   * @param i Index of a walker
   * @return Its intended tempering threshold in kcal/mol.
   */
  private double getTemperingThreshold(int i) {
    return group.temperingThreshold.length > 1 ? group.temperingThreshold[i]
        : group.temperingThreshold[0];
  }

  /**
   * Returns the tempering parameter associated with a walker.
   *
   * @param i Index of a walker
   * @return Its intended tempering parameter in kBT.
   */
  private double getTemperingParameter(int i) {
    return group.temperingRate.length > 1 ? group.temperingRate[i] : group.temperingRate[0];
  }

  /**
   * Sets the number of time steps between OST counts.
   *
   * @return Returns the interval between OST counts.
   */
  public int getCountInterval() {
    return group.countInterval;
  }

  public void setCountInterval(int countInterval) {
    group.countInterval = countInterval;
  }

  /**
   * Sets the initial Gaussian bias magnitude in kcal/mol.
   *
   * @return Returns the Bias magnitude.
   */
  public double[] getBiasMag() {
    return group.biasMag;
  }

  public void setBiasMag(double[] biasMag) {
    group.biasMag = biasMag;
  }

  /**
   * Enforces that each walker maintains their own histogram.
   *
   * @return Returns true if each Walker has their own histogram.
   */
  public boolean isIndependentWalkers() {
    return group.independentWalkers;
  }

  public void setIndependentWalkers(boolean independentWalkers) {
    group.independentWalkers = independentWalkers;
  }

  /**
   * The Dama et al tempering rate parameter, in multiples of kBT.
   *
   * @return Returns the tempering rate.
   */
  public double[] getTemperingRate() {
    return group.temperingRate;
  }

  public void setTemperingRate(double[] temperingRate) {
    group.temperingRate = temperingRate;
  }

  /**
   * The tempering threshold/offset in kcal/mol.
   *
   * @return Returns the tempering threshold.
   */
  public double[] getTemperingThreshold() {
    return group.temperingThreshold;
  }

  public void setTemperingThreshold(double[] temperingThreshold) {
    group.temperingThreshold = temperingThreshold;
  }

  /**
   * The Monte Carlo scheme can use a hard wall that rejects any sample (Lambda, dU/dL) located in an
   * empty histogram bin.
   *
   * @return Returns true if the MC-OST hard wall constraint is set.
   */
  public boolean isMcHardWall() {
    return mcGroup.mcHardWall;
  }

  public void setMcHardWall(boolean mcHardWall) {
    mcGroup.mcHardWall = mcHardWall;
  }

  /**
   * The number of steps to take for each MD trajectory for MC-OST.
   *
   * @return Returns the number of MD steps for each MC-OST round.
   */
  public int getMcMDSteps() {
    return mcGroup.mcMDSteps;
  }

  public void setMcMDSteps(int mcMDSteps) {
    mcGroup.mcMDSteps = mcMDSteps;
  }

  /**
   * The standard deviation for lambda moves.
   *
   * @return Returns the MC lambda trial move standard deviations.
   */
  public double getMcLambdaStdDev() {
    return mcGroup.mcLambdaStdDev;
  }

  public void setMcLambdaStdDev(double mcLambdaStdDev) {
    mcGroup.mcLambdaStdDev = mcLambdaStdDev;
  }

  /**
   * Only write out snapshots if lambda is greater than the value specified.
   *
   * @return Returns the lambda write-out threshold.
   */
  public double getLambdaWriteOut() {
    return group.lambdaWriteOut;
  }

  public void setLambdaWriteOut(double lambdaWriteOut) {
    group.lambdaWriteOut = lambdaWriteOut;
  }

  /**
   * Collection of Orthogonal Space Tempering Options.
   */
  private static class OSTOptionGroup {

    /** -c or --count Sets the number of time steps between OST counts. */
    @Option(
        names = {"-C", "--count"},
        paramLabel = "10",
        defaultValue = "10",
        description = "Time steps between MD Orthogonal Space counts.")
    private int countInterval;

    /** --bM or --biasMag sets the initial Gaussian bias magnitude in kcal/mol. */
    @Option(
        names = {"--bM", "--biasMag"},
        paramLabel = "0.05",
        defaultValue = "0.05",
        split = ",",
        description =
            "Orthogonal Space Gaussian bias magnitude (kcal/mol); RepEx OST uses a comma-separated list.")
    private double[] biasMag;

    /** --iW or --independentWalkers enforces that each walker maintains their own histogram. */
    @Option(
        names = {"--iW", "--independentWalkers"},
        defaultValue = "false",
        description = "Enforces that each walker maintains their own histogram. ")
    private boolean independentWalkers;

    /** --tp or --temperingRate sets the Dama et al tempering rate parameter, in multiples of kBT. */
    @Option(
        names = {"--tp", "--temperingRate"},
        paramLabel = "4.0",
        defaultValue = "4.0",
        split = ",",
        description =
            "Tempering rate parameter in multiples of kBT; RepEx OST uses a comma-separated list.")
    private double[] temperingRate;

    /** --tth or --temperingThreshold sets the tempering threshold/offset in kcal/mol. */
    @Option(
        names = {"--tth", "--temperingThreshold"},
        paramLabel = "20*bias",
        defaultValue = "-1",
        split = ",",
        description = "Tempering threshold in kcal/mol; RepEx OST uses a comma-separated list.")
    private double[] temperingThreshold;

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
  }

  /**
   * Collection of Monte Carlo Orthogonal Space Tempering Options.
   */
  private static class MCOSTOptionGroup {

    /** --mc or --monteCarlo sets the Monte Carlo scheme for Orthogonal Space Tempering. */
    @Option(
        names = {"--mc", "--monteCarlo"},
        defaultValue = "false",
        description = "Specify use of Monte Carlo OST")
    private boolean monteCarlo;

    /**
     * --mcHW or --mcHardWall sets the Monte Carlo scheme to use a hard wall that rejects any sample
     * (Lambda, dU/dL) located in an empty histogram bin.
     */
    @Option(
        names = {"--mcHW", "--mcHardWall"},
        defaultValue = "false",
        description = "Monte Carlo OST hard wall constraint.")
    private boolean mcHardWall;

    /** --mcmD or --mcMDSteps Sets the number of steps to take for each MD trajectory for MC-OST. */
    @Option(
        names = {"--mcMD", "--mcMDSteps"},
        paramLabel = "100",
        defaultValue = "100",
        description = "Number of dynamics steps to take for each MD trajectory for Monte Carlo OST")
    private int mcMDSteps;

    /** --mcL or --mcLambdaStdDev Sets the standard deviation for lambda moves. */
    @Option(
        names = {"--mcL", "--mcLambdaStdDev"},
        paramLabel = "0.01",
        defaultValue = "0.01",
        description = "Standard deviation for lambda move.")
    private double mcLambdaStdDev;

    /** --ts or --twoStep MC Orthogonal Space sampling using separate lambda and MD moves. */
    @Option(
        names = {"--ts", "--twoStep"},
        defaultValue = "false",
        description = "MC Orthogonal Space sampling using separate lambda and MD moves.")
    private boolean twoStep;
  }

}
