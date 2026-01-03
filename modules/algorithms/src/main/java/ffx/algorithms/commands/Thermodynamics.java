//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import edu.rit.pj.Comm;
import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.cli.BarostatOptions;
import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.cli.LambdaParticleOptions;
import ffx.algorithms.cli.MultiDynamicsOptions;
import ffx.algorithms.cli.OSTOptions;
import ffx.algorithms.cli.RandomUnitCellOptions;
import ffx.algorithms.cli.ThermodynamicsOptions;
import ffx.algorithms.cli.ThermodynamicsOptions.ThermodynamicsAlgorithm;
import ffx.algorithms.thermodynamics.MonteCarloOST;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering;
import ffx.crystal.CrystalPotential;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.cli.AlchemicalOptions;
import ffx.potential.cli.TopologyOptions;
import ffx.potential.cli.WriteoutOptions;
import ffx.utilities.FFXBinding;
import org.apache.commons.configuration2.Configuration;
import org.apache.commons.io.FilenameUtils;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static java.lang.String.format;

/**
 * The Thermodynamics script uses the Transition-Tempered Orthogonal Space Random Walk
 * algorithm to estimate a free energy.
 * <br>
 * Usage:
 * <br>
 * ffxc Thermodynamics [options] &lt;filename&gt; [file2...]
 */
@Command(description = " Use the Transition-Tempered Orthogonal Space Random Walk algorithm to estimate a free energy.", name = "Thermodynamics")
public class Thermodynamics extends AlgorithmsCommand {

  @Mixin
  public DynamicsOptions dynamicsOptions;

  @Mixin
  public BarostatOptions barostatOptions;

  @Mixin
  public RandomUnitCellOptions randomSymopOptions;

  @Mixin
  public AlchemicalOptions alchemicalOptions;

  @Mixin
  public TopologyOptions topologyOptions;

  @Mixin
  public WriteoutOptions writeoutOptions;

  @Mixin
  public ThermodynamicsOptions thermodynamicsOptions;

  @Mixin
  public OSTOptions ostOptions;

  @Mixin
  public LambdaParticleOptions lambdaParticleOptions;

  @Mixin
  public MultiDynamicsOptions multiDynamicsOptions;

  /**
   * -v or --verbose  Log additional information (primarily for MC-OST).
   */
  @Option(names = {"-v", "--verbose"},
      description = "Log additional information (primarily for MC-OST).")
  public boolean verbose = false;

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "The atomic coordinate file in PDB or XYZ format.")
  public List<String> filenames = null;

  public MolecularAssembly[] topologies;
  public CrystalPotential potential;
  public OrthogonalSpaceTempering orthogonalSpaceTempering = null;
  public Configuration additionalProperties;

  /**
   * Sets an optional Configuration with additional properties.
   *
   * @param additionalProps Additional properties configuration
   */
  public void setProperties(Configuration additionalProps) {
    this.additionalProperties = additionalProps;
  }

  /**
   * Thermodynamics Constructor.
   */
  public Thermodynamics() {
    super();
  }

  /**
   * Thermodynamics Constructor.
   *
   * @param binding The Binding to use.
   */
  public Thermodynamics(FFXBinding binding) {
    super(binding);
  }

  /**
   * Thermodynamics constructor that sets the command line arguments.
   *
   * @param args Command line arguments.
   */
  public Thermodynamics(String[] args) {
    super(args);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public Thermodynamics run() {

    // Begin boilerplate "make a topology" code.
    if (!init()) {
      return this;
    }

    // Determine the number of topologies to be read and allocate the array.
    int numTopologies = topologyOptions.getNumberOfTopologies(filenames);
    int threadsPerTopology = topologyOptions.getThreadsPerTopology(numTopologies);
    topologies = new MolecularAssembly[numTopologies];

    // Turn on computation of lambda derivatives if softcore atoms exist.
    alchemicalOptions.setAlchemicalProperties();
    topologyOptions.setAlchemicalProperties(numTopologies);

    Comm world = Comm.world();
    int size = world.size();
    int rank = (size > 1) ? world.rank() : 0;

    // Segment of code for MultiDynamics and OST.
    List<File> structureFiles = new ArrayList<>();
    for (String filename : filenames) {
      File file = new File(FilenameUtils.normalize(filename));
      structureFiles.add(file);
    }

    File firstStructure = structureFiles.get(0);
    String filePathNoExtension = firstStructure.getAbsolutePath().replaceFirst("\\.[^.]+$", "");
    File histogramRestart = new File(filePathNoExtension + ".his");

    // For a multi-process job, try to get the restart files from rank sub-directories.
    String withRankName = filePathNoExtension;

    if (size > 1) {
      List<File> rankedFiles = new ArrayList<>(numTopologies);
      String rankDirName = FilenameUtils.getFullPath(filePathNoExtension);
      rankDirName = format("%s%d", rankDirName, rank + multiDynamicsOptions.getFirstDir());
      File rankDirectory = new File(rankDirName);
      if (!rankDirectory.exists()) {
        rankDirectory.mkdir();
      }
      rankDirName = rankDirName + File.separator;
      withRankName = format("%s%s", rankDirName, FilenameUtils.getName(filePathNoExtension));

      for (File structureFile : structureFiles) {
        rankedFiles.add(new File(format("%s%s", rankDirName,
            FilenameUtils.getName(structureFile.getName()))));
      }
      structureFiles = rankedFiles;
    }

    File lambdaRestart = new File(withRankName + ".lam");
    File dyn = new File(withRankName + ".dyn");
    if (ostOptions.getIndependentWalkers()) {
      histogramRestart = new File(withRankName + ".his");
    }

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
        topologies[i] = multiDynamicsOptions.openFile(algorithmFunctions, topologyOptions,
            threadsPerTopology, filenames.get(i), i, alchemicalOptions, structureFiles.get(i), rank);
      }
    }

    StringBuilder sb = new StringBuilder("\n Running ");

    ThermodynamicsAlgorithm algorithm = thermodynamicsOptions.getAlgorithm();
    double initLambda = alchemicalOptions.getInitialLambda(size, rank, true);
    if (algorithm == ThermodynamicsAlgorithm.OST) {
      sb.append("Orthogonal Space Tempering");
    } else if (algorithm == ThermodynamicsAlgorithm.FIXED) {
      sb.append("Fixed Lambda Sampling at Window L=").append(format("%5.3f ", initLambda));
    } else if (algorithm == ThermodynamicsAlgorithm.NEQ) {
      sb.append("Non-Equilibrium Sampling");
    } else {
      logger.severe(" Unknown Thermodynamics Algorithm " + algorithm);
    }
    sb.append(" for ");

    potential = (CrystalPotential) topologyOptions.assemblePotential(topologies, sb);
    logger.info(sb.toString());

    LambdaInterface lambdaInterface = (LambdaInterface) potential;
    boolean lamExists = lambdaRestart.exists();
    double[] x = new double[potential.getNumberOfVariables()];
    potential.getCoordinates(x);
    lambdaInterface.setLambda(initLambda);
    potential.energy(x, true);

    if (numTopologies == 1) {
      randomSymopOptions.randomize(topologies[0]);
    }

    multiDynamicsOptions.distribute(topologies, potential, algorithmFunctions, rank, size);

    if (algorithm == ThermodynamicsAlgorithm.OST) {
      orthogonalSpaceTempering = ostOptions.constructOST(potential, lambdaRestart, histogramRestart, topologies[0],
              additionalProperties, dynamicsOptions, thermodynamicsOptions, lambdaParticleOptions,
              algorithmListener, !multiDynamicsOptions.isSynchronous());
      if (!lamExists) {
        orthogonalSpaceTempering.setLambda(initLambda);
      }
      // Can be either the OST or a Barostat on top of it.
      CrystalPotential ostPotential = ostOptions.applyAllOSTOptions(orthogonalSpaceTempering, topologies[0],
          dynamicsOptions, barostatOptions);
      if (ostOptions.isMonteCarlo()) {
        MonteCarloOST mcOST = ostOptions.setupMCOST(orthogonalSpaceTempering, topologies, ostPotential,
            dynamicsOptions, thermodynamicsOptions, verbose, dyn, algorithmListener);
        ostOptions.beginMCOST(mcOST, dynamicsOptions, thermodynamicsOptions);
      } else {
        ostOptions.beginMDOST(orthogonalSpaceTempering, topologies, ostPotential, dynamicsOptions,
            writeoutOptions, thermodynamicsOptions, dyn, algorithmListener);
      }
      logger.info(" Done running OST sampling.");
    } else if (algorithm == ThermodynamicsAlgorithm.FIXED) {
      orthogonalSpaceTempering = null;
      potential = barostatOptions.checkNPT(topologies[0], potential);
      thermodynamicsOptions.runFixedAlchemy(topologies, potential, dynamicsOptions, writeoutOptions, dyn, algorithmListener);
      logger.info(" Done running fixed lambda sampling.");
    } else if (algorithm == ThermodynamicsAlgorithm.NEQ) {
      orthogonalSpaceTempering = null;
      potential = barostatOptions.checkNPT(topologies[0], potential);
      thermodynamicsOptions.runNEQ(topologies, potential, dynamicsOptions, writeoutOptions, dyn, algorithmListener);
      logger.info(" Done running non-equilibrium sampling.");
    } else {
      logger.severe(" Unknown Thermodynamics Algorithm " + algorithm);
    }

    return this;
  }

  public OrthogonalSpaceTempering getOST() {
    return orthogonalSpaceTempering;
  }

  public CrystalPotential getPotential() {
    return potential;
  }

  @Override
  public List<Potential> getPotentials() {
    List<Potential> potentials;
    if (orthogonalSpaceTempering == null) {
      if (potential == null) {
        potentials = Collections.emptyList();
      } else {
        potentials = Collections.singletonList((Potential) potential);
      }
    } else {
      potentials = Collections.singletonList((Potential) orthogonalSpaceTempering);
    }
    return potentials;
  }

  @Override
  public boolean destroyPotentials() {
    return getPotentials().stream().allMatch(potential -> potential.destroy());
  }
}
