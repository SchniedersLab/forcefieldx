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

import static java.lang.String.format;

import ffx.algorithms.AlgorithmFunctions;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.optimize.RotamerOptimization;
import ffx.crystal.CrystalPotential;
import ffx.numerics.Potential;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.QuadTopologyEnergy;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.cli.AlchemicalOptions;
import ffx.potential.cli.TopologyOptions;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that can create multiple walkers, such as
 * multi-walker OST. Should be kept agnostic to whether it is an MD-based algorithm, or some other
 * flavor of Monte Carlo.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class MultiDynamicsOptions {

  private static final Logger logger = Logger.getLogger(MultiDynamicsOptions.class.getName());

  /** -y or --synchronous sets synchronous walker communication (not recommended) */
  @Option(
      names = {"-y", "--synchronous"},
      defaultValue = "false",
      description = "Walker communication is synchronous")
  private boolean synchronous;

  /**
   * -dw or --distributeWalkers allows walkers to start from multiple conformations; AUTO picks up
   * per-walker conformations as filename.pdb_(walker number), and specifying a residue starts a
   * rotamer optimization to generate side-chain configurations to start from.
   */
  @Option(
      names = {"--dw", "--distributeWalkers"},
      paramLabel = "OFF",
      defaultValue = "OFF",
      description =
          "AUTO: Pick up per-walker configurations as [filename.pdb]_[num], or specify a residue to distribute on.")
  private String distributeWalkersString;

  /**
   * If residues selected for distributing initial configurations, performs many-body optimization
   * for this distribution.
   *
   * @param molecularAssemblies an array of {@link ffx.potential.MolecularAssembly} objects.
   * @param crystalPotential Overall CrystalPotential in use.
   * @param algorithmFunctions a {@link ffx.algorithms.AlgorithmFunctions} object.
   * @param rank a int.
   * @param worldSize a int.
   */
  public void distribute(
      MolecularAssembly[] molecularAssemblies,
      CrystalPotential crystalPotential,
      AlgorithmFunctions algorithmFunctions,
      int rank,
      int worldSize) {
    int ntops = molecularAssemblies.length;
    Potential[] energies = new Potential[ntops];
    for (int i = 0; i < ntops; i++) {
      energies[i] = molecularAssemblies[i].getPotentialEnergy();
    }
    distribute(
        molecularAssemblies, energies, crystalPotential, algorithmFunctions, rank, worldSize);
  }

  /**
   * If residues selected for distributing initial configurations, performs many-body optimization
   * for this distribution.
   *
   * @param molecularAssemblies an array of {@link ffx.potential.MolecularAssembly} objects.
   * @param potentials ForceFieldEnergy for each topology.
   * @param crystalPotential Overall CrystalPotential in use.
   * @param algorithmFunctions a {@link ffx.algorithms.AlgorithmFunctions} object.
   * @param rank a int.
   * @param worldSize a int.
   */
  public void distribute(
      MolecularAssembly[] molecularAssemblies,
      Potential[] potentials,
      CrystalPotential crystalPotential,
      AlgorithmFunctions algorithmFunctions,
      int rank,
      int worldSize) {
    if (!distributeWalkersString.equalsIgnoreCase("AUTO")
        && !distributeWalkersString.equalsIgnoreCase("OFF")) {
      logger.info(" Distributing walker conformations.");
      int nSys = molecularAssemblies.length;
      assert nSys == potentials.length;
      switch (nSys) {
        case 1:
          optStructure(
              molecularAssemblies[0], crystalPotential, algorithmFunctions, rank, worldSize);
          break;
        case 2:
          DualTopologyEnergy dte = (DualTopologyEnergy) crystalPotential;
          if (dte.getNumSharedVariables() == dte.getNumberOfVariables()) {
            logger.info(" Generating starting structures based on dual-topology:");
            optStructure(molecularAssemblies[0], dte, algorithmFunctions, rank, worldSize);
          } else {
            logger.info(
                " Generating separate starting structures for each topology of the dual toplogy:");
            optStructure(
                molecularAssemblies[0], potentials[0], algorithmFunctions, rank, worldSize);
            optStructure(
                molecularAssemblies[1], potentials[1], algorithmFunctions, rank, worldSize);
          }
          break;
        case 4:
          QuadTopologyEnergy qte = (QuadTopologyEnergy) crystalPotential;
          optStructure(
              molecularAssemblies[0], qte.getDualTopA(), algorithmFunctions, rank, worldSize);
          optStructure(
              molecularAssemblies[3], qte.getDualTopB(), algorithmFunctions, rank, worldSize);
          break;
        // Oct-topology is deprecated on account of not working as intended.
        default:
          logger.severe(" First: must have 1, 2, or 4 topologies.");
          break;
      }
    } else {
      logger.finer(" Skipping RO-based distribution of initial configurations.");
    }
  }

  /**
   * Synchronous walker communication.
   *
   * <p>isSynchronous.
   *
   * @return a boolean.
   */
  public boolean isSynchronous() {
    return synchronous;
  }

  public void setSynchronous(boolean synchronous) {
    this.synchronous = synchronous;
  }

  /**
   * Opens a file and processes it. Extends the behavior of AlchemicalOptions.openFile by permitting
   * use of a rank-dependent File.
   *
   * @param algorithmFunctions AlgorithmFunctions object.
   * @param topologyOptions Topology Options.
   * @param threadsPer Threads to use per system.
   * @param toOpen Filename to open.
   * @param topNum Number of the topology to open.
   * @param alchemicalOptions Alchemical Options.
   * @param rank Rank in the world communicator.
   * @param structureFile a {@link java.io.File} object.
   * @return a {@link ffx.potential.MolecularAssembly} object.
   */
  public MolecularAssembly openFile(
      AlgorithmFunctions algorithmFunctions,
      TopologyOptions topologyOptions,
      int threadsPer,
      String toOpen,
      int topNum,
      AlchemicalOptions alchemicalOptions,
      File structureFile,
      int rank) {
    boolean autoDist = distributeWalkersString.equalsIgnoreCase("AUTO");

    if (autoDist) {
      String openName = format("%s_%d", toOpen, rank + 1);
      File testFile = new File(openName);
      if (testFile.exists()) {
        toOpen = openName;
      } else {
        logger.warning(format(" File %s does not exist; using default %s", openName, toOpen));
      }
    }
    MolecularAssembly assembly =
        alchemicalOptions.openFile(algorithmFunctions, topologyOptions, threadsPer, toOpen, topNum);
    assembly.setFile(structureFile);
    return assembly;
  }

  /**
   * Parses --dw into optimization tokens if it's not "OFF", "AUTO", or null.
   *
   * @return An array of Strings from splitting the distributed flag.
   */
  private String[] parseDistributed() {
    if (distributeWalkersString.equalsIgnoreCase("OFF")
        || distributeWalkersString.equalsIgnoreCase("AUTO")
        || distributeWalkersString.isEmpty()) {
      return null;
    }
    return distributeWalkersString.split("\\.");
  }

  /**
   * Distribute side-chain conformations of molecularAssembly.
   *
   * @param molecularAssembly To distribute
   * @param potential Potential to use
   */
  private void optStructure(
      MolecularAssembly molecularAssembly,
      Potential potential,
      AlgorithmFunctions algorithmFunctions,
      int rank,
      int worldSize) {
    RotamerLibrary rLib = new RotamerLibrary(false);
    String[] distribRes = parseDistributed();

    if (distribRes == null || distribRes.length == 0) {
      throw new IllegalArgumentException(
          " Programming error: Must have list of residues to split on!");
    }

    LambdaInterface lambdaInterface =
        (potential instanceof LambdaInterface) ? (LambdaInterface) potential : null;
    double initLam = -1.0;
    if (lambdaInterface != null) {
      initLam = lambdaInterface.getLambda();
      lambdaInterface.setLambda(0.5);
    }

    Pattern chainMatcher = Pattern.compile("^([a-zA-Z])?([0-9]+)$");

    List<Residue> residueList = new ArrayList<>(distribRes.length);

    for (String ts : distribRes) {
      Matcher m = chainMatcher.matcher(ts);
      Character chainID;
      int resNum;
      if (m.find()) {
        if (m.groupCount() == 2) {
          chainID = m.group(1).charAt(0);
          resNum = Integer.parseInt(m.group(2));
        } else {
          chainID = ' ';
          resNum = Integer.parseInt(m.group(1));
        }
      } else {
        logger.warning(format(" Could not parse %s as a valid residue!", ts));
        continue;
      }
      logger.info(format(" Looking for chain %c residue %d", chainID, resNum));

      for (Polymer p : molecularAssembly.getChains()) {
        if (p.getChainID() == chainID) {
          for (Residue r : p.getResidues()) {
            if (r.getResidueNumber() == resNum && r.getRotamers(rLib) != null) {
              residueList.add(r);
            }
          }
        }
      }
    }

    if (residueList.isEmpty()) {
      throw new IllegalArgumentException(" No valid entries for distWalkers!");
    }

    AlgorithmListener alist = algorithmFunctions.getDefaultListener();
    RotamerOptimization ropt = new RotamerOptimization(molecularAssembly, potential, alist);
    ropt.setRotamerLibrary(rLib);

    ropt.setThreeBodyEnergy(false);

    CompositeConfiguration properties = molecularAssembly.getProperties();
    if (!properties.containsKey("ro-ensembleNumber")
        && !properties.containsKey("ro-ensembleEnergy")) {
      logger.info(format(" Setting ensemble to default of number of walkers %d", worldSize));
      ropt.setEnsemble(worldSize);
    }

    ropt.setPrintFiles(false);
    ropt.setResiduesIgnoreNull(residueList);

    RotamerLibrary.measureRotamers(residueList, false);

    boolean lazyMat = properties.getBoolean("ro-lazyMatrix", false);
    if (!lazyMat) {
      properties.clearProperty("ro-lazyMatrix");
      properties.addProperty("ro-lazyMatrix", true);
    }

    ropt.optimize(RotamerOptimization.Algorithm.ALL);
    ropt.setCoordinatesToEnsemble(rank);

    // One final energy call to ensure the coordinates are properly set at the
    // end of rotamer optimization.
    double[] xyz = new double[potential.getNumberOfVariables()];
    potential.getCoordinates(xyz);
    logger.info(" Final Optimized Energy:");
    potential.energy(xyz, true);

    if (lambdaInterface != null) {
      lambdaInterface.setLambda(initLam);
    }

    if (!lazyMat) {
      properties.clearProperty("ro-lazyMatrix");
    }
  }

  /**
   * Allows walkers to start from multiple conformations; AUTO picks up per-walker conformations as
   * filename.pdb_(walker number), and specifying a residue starts a rotamer optimization to
   * generate side-chain configurations to start from.
   */
  public String getDistributeWalkersString() {
    return distributeWalkersString;
  }

  public void setDistributeWalkersString(String distributeWalkersString) {
    this.distributeWalkersString = distributeWalkersString;
  }
}
