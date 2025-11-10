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
package ffx.xray.commands.test;

import edu.rit.pj.Comm;
import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.cli.LambdaParticleOptions;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.algorithms.thermodynamics.HistogramData;
import ffx.algorithms.thermodynamics.LambdaData;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering.OptimizationParameters;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
import ffx.realspace.cli.RealSpaceOptions;
import ffx.utilities.FFXBinding;
import ffx.xray.DiffractionData;
import ffx.xray.RefinementEnergy;
import ffx.xray.refine.RefinementMode;
import ffx.xray.cli.XrayOptions;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * The Alchemical Changes script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.Alchemical [options] &lt;filename&gt;
 */
@Command(description = " Simulated annealing on an X-ray target.", name = "xray.test.Alchemical")
public class Alchemical extends AlgorithmsCommand {

  @Mixin
  private DynamicsOptions dynamicsOptions;

  @Mixin
  private LambdaParticleOptions lambdaParticleOptions;

  @Mixin
  private XrayOptions xrayOptions;

  private static final Logger logger = Logger.getLogger(RealSpaceOptions.class.getName());

  /**
   * -I or --onlyIons sets whether or not ion positions are optimized (default is false; must set at least one of either '-W' or '-I') (only one type of ion is chosen).
   */
  @Option(names = {"-I", "--onlyIons"}, paramLabel = "false",
      description = "Set to only optimize ions (of a single type).")
  private boolean onlyIons = false;

  /**
   * --itype or --iontype Specify which ion to run optimization on. If none is specified, default behavior chooses the first ion found in the PDB file.
   */
  @Option(names = {"--itype", "--iontype"}, paramLabel = "null",
      description = "Specify which ion to run optimization on. If none is specified, default behavior chooses the first ion found in the PDB file.")
  private String[] ionType = null;

  /**
   * -N or --neutralize Adds more of the selected ion in order to neutralize the crystal's charge.
   */
  @Option(names = {"-N", "--neutralize"}, paramLabel = "false",
      description = "Neutralize the crystal's charge by adding more of the selected ion")
  private boolean neutralize = false;

  /**
   * -W or --onlyWater sets whether or not water positions are optimized (default is false; must set at least one of either '-W' or '-I').
   */
  @Option(names = {"-W", "--onlyWater"}, paramLabel = "false",
      description = "Set to only optimize water.")
  private boolean onlyWater = false;

  /**
   * The refinement mode to use.
   */
  private RefinementMode refinementMode = RefinementMode.COORDINATES;

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
  private List<String> filenames;

  // Value of Lambda.
  private double lambda = 1.0;

  // ThermostatEnum [ ADIABATIC, BERENDSEN, BUSSI ]
  private ThermostatEnum thermostat = ThermostatEnum.ADIABATIC;

  // IntegratorEnum [ BEEMAN, RESPA, STOCHASTIC ]
  private IntegratorEnum integrator = IntegratorEnum.STOCHASTIC;

  // File type of coordinate snapshots to write out.
  private String fileType = "PDB";

  // OST
  private boolean runOST = true;

  // Reset velocities (ignored if a restart file is given)
  private boolean initVelocities = true;

  private OrthogonalSpaceTempering orthogonalSpaceTempering;

  /**
   * Alchemical constructor.
   */
  public Alchemical() {
    super();
  }

  /**
   * Alchemical constructor that sets the command line arguments.
   * @param args Command line arguments.
   */
  public Alchemical(String[] args) {
    super(args);
  }

  /**
   * Alchemical constructor.
   * @param binding The Binding to use.
   */
  public Alchemical(FFXBinding binding) {
    super(binding);
  }

  @Override
  public Alchemical run() {

    if (!init()) {
      return this;
    }
    dynamicsOptions.init();
    xrayOptions.init();
    System.setProperty("lambdaterm", "true");

    String modelfilename;
    MolecularAssembly[] assemblies;
    if (filenames != null && filenames.size() > 0) {
      assemblies = algorithmFunctions.openAll(filenames.get(0));
      activeAssembly = assemblies[0];
      modelfilename = filenames.get(0);
    } else if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    } else {
      modelfilename = activeAssembly.getFile().getAbsolutePath();
      assemblies = new MolecularAssembly[]{activeAssembly};
    }

    logger.info("\n Running Alchemical Changes on " + modelfilename);

    File structureFile = new File(FilenameUtils.normalize(modelfilename));
    structureFile = new File(structureFile.getAbsolutePath());
    String baseFilename = FilenameUtils.removeExtension(structureFile.getName());
    File histogramRestart = new File(baseFilename + ".his");
    File lambdaRestart = new File(baseFilename + ".lam");
    File dyn = new File(baseFilename + ".dyn");

    Comm world = Comm.world();
    int size = world.size();
    int rank = 0;
    double[] energyArray = new double[world.size()];
    for (int i = 0; i < world.size(); i++) {
      energyArray[i] = Double.MAX_VALUE;
    }

    // For a multi-process job, try to get the restart files from rank sub-directories.
    if (size > 1) {
      rank = world.rank();
      File rankDirectory = new File(
          structureFile.getParent() + File.separator + Integer.toString(rank));
      if (!rankDirectory.exists()) {
        rankDirectory.mkdir();
      }
      lambdaRestart = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam");
      dyn = new File(rankDirectory.getPath() + File.separator + baseFilename + ".dyn");
      structureFile = new File(rankDirectory.getPath() + File.separator + structureFile.getName());
    }
    if (!dyn.exists()) {
      dyn = null;
    }

    // Set built atoms active/use flags to true (false for other atoms).
    Atom[] atoms = activeAssembly.getAtomArray();

    // Get a reference to the first system's ForceFieldEnergy.
    ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy();
    forceFieldEnergy.setPrintOnFailure(false, false);

    // Configure all atoms to:
    // 1) be used in the potential
    // 2) be inactive (i.e. cannot move)
    // 3) not be controlled by the lambda state variable.
    for (int i = 0; i <= atoms.length; i++) {
      Atom ai = atoms[i - 1];
      ai.setUse(true);
      ai.setActive(false);
      ai.setApplyLambda(false);
    }

    double crystalCharge = activeAssembly.getCharge(true);
    logger.info(" Overall crystal charge: " + crystalCharge);
    List<MSNode> ions = assemblies[0].getIons();
    List<MSNode> water = assemblies[0].getWater();

    // Consider the option of creating a composite lambda gradient from vapor phase to crystal phase
    if (!onlyWater) {
      logger.info("Doing ions.");
      if (ions == null || ions.size() == 0) {
        logger.info("\n Please add an ion to the PDB file to scan with.");
        return this;
      }
      for (MSNode msNode : ions) {
        // Check if ionType is specified and matches this ion's atom names
        boolean ionSelected = false;
        if (ionType != null) {
          logger.info("Selecting ion.");
          List<Atom> atomList = msNode.getAtomList();
          if (!atomList.isEmpty()) {
            String atomName = atomList.get(0).getName();
            for (String targetIonType : ionType) {
              if (atomName.equals(targetIonType)) {
                ionSelected = true;
                break;
              }
            }
          }
        }

        if (ionSelected) {
          logger.info("Ion has been selected.");
          for (Atom atom : msNode.getAtomList()) {
            System.out.println("Activating ions");
            atom.setUse(true);
            atom.setActive(true);
            atom.setApplyLambda(true);
            logger.info(" Alchemical atom: " + atom.toString());
          }
        } else {
          logger.info("Ion has not been selected.");
          if (neutralize) {
            logger.info("Neutralizing crystal.");
            double ionCharge = 0;
            for (Atom atom : msNode.getAtomList()) {
              ionCharge += atom.getMultipoleType().getCharge();
            }
            logger.info("Ion charge is: " + Double.toString(ionCharge));
            int numIons = (int) (-1 * (Math.ceil(crystalCharge / ionCharge)));
            if (numIons > 0) {
              List<Atom> atomList = msNode.getAtomList();
              String atomName = atomList.isEmpty() ? "" : atomList.get(0).getName();
              logger.info(numIons + " " + atomName
                  + " ions needed to neutralize the crystal.");
              ionType = new String[]{atomName};
              for (Atom atom : msNode.getAtomList()) {
                atom.setUse(true);
                atom.setActive(true);
                atom.setApplyLambda(true);
                logger.info(" Alchemical atom: " + atom.toString());
              }
            }
          } else {
            List<Atom> atomList = msNode.getAtomList();
            String atomName = atomList.isEmpty() ? "" : atomList.get(0).getName();
            ionType = new String[]{atomName};
            for (Atom atom : msNode.getAtomList()) {
              atom.setUse(true);
              atom.setActive(true);
              atom.setApplyLambda(true);
              logger.info(" Alchemical atom: " + atom.toString());
            }
          }
        }
      }
    }

    // Lambdize water for position optimization, if this option was set to true
    if (onlyWater) {
      for (MSNode msNode : water) {
        for (Atom atom : msNode.getAtomList()) {
          // Scan with the last ion in the file.
          atom.setUse(true);
          atom.setActive(true);
          atom.setApplyLambda(true);
          logger.info(" Water atom:      " + atom.toString());
        }
      }
    }

    // Load parsed X-ray properties.
    CompositeConfiguration properties = assemblies[0].getProperties();
    xrayOptions.setProperties(parseResult, properties);

    DiffractionData diffractionData = xrayOptions.getDiffractionData(filenames, assemblies, properties);
    RefinementEnergy refinementEnergy = xrayOptions.toXrayEnergy(diffractionData);

    double[] x = new double[refinementEnergy.getNumberOfVariables()];
    x = refinementEnergy.getCoordinates(x);
    refinementEnergy.energy(x, true);

    refinementEnergy.setLambda(lambda);

    CompositeConfiguration props = assemblies[0].getProperties();

    HistogramData histogramData = HistogramData.readHistogram(histogramRestart);
    if (lambdaRestart == null) {
      String filename = histogramRestart.toString().replaceFirst("\\.his$", ".lam");
      lambdaRestart = new File(filename);
    }
    LambdaData lambdaData = LambdaData.readLambdaData(lambdaRestart);

    OrthogonalSpaceTempering orthogonalSpaceTempering = new OrthogonalSpaceTempering(refinementEnergy,
        refinementEnergy, histogramData, lambdaData, props, dynamicsOptions, lambdaParticleOptions, algorithmListener);

    orthogonalSpaceTempering.setLambda(lambda);
    orthogonalSpaceTempering.getOptimizationParameters().setOptimization(true, activeAssembly);
    // Create the MolecularDynamics instance.
    MolecularDynamics molDyn = new MolecularDynamics(assemblies[0],
        orthogonalSpaceTempering, algorithmListener, thermostat, integrator);

    algorithmFunctions.energy(assemblies[0]);

    molDyn.dynamic(dynamicsOptions.getSteps(), dynamicsOptions.getDt(), dynamicsOptions.getReport(),
        dynamicsOptions.getWrite(), dynamicsOptions.getTemperature(), true,
        fileType, dynamicsOptions.getWrite(), dyn);

    logger.info(" Searching for low energy coordinates");
    OptimizationParameters opt = orthogonalSpaceTempering.getOptimizationParameters();
    double[] lowEnergyCoordinates = opt.getOptimumCoordinates();
    double currentOSTOptimum = opt.getOptimumEnergy();
    if (lowEnergyCoordinates != null) {
      forceFieldEnergy.setCoordinates(lowEnergyCoordinates);
    } else {
      logger.info(" OST stage did not succeed in finding a minimum.");
    }

    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    return orthogonalSpaceTempering == null ? Collections.emptyList() :
        Collections.singletonList(orthogonalSpaceTempering);
  }
}
