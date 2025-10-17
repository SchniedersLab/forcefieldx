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

import edu.rit.pj.Comm;
import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.cli.RepExOptions;
import ffx.algorithms.dynamics.MDEngine;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.dynamics.MolecularDynamicsOpenMM;
import ffx.algorithms.dynamics.PhReplicaExchange;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.cli.AtomSelectionOptions;
import ffx.potential.cli.WriteoutOptions;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XPHFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.FFXBinding;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import picocli.CommandLine.Command;
import picocli.CommandLine.Mixin;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static java.lang.String.format;

/**
 * The PhDynamics command implements constant pH molecular dynamics.
 * <br>
 * Usage:
 * <br>
 * ffxc PhDynamics [options] &lt;filename&gt; [file2...]
 */
@Command(description = " Run constant pH dynamics on a system.", name = "PhDynamics")
public class PhDynamics extends AlgorithmsCommand {

  @Mixin
  private AtomSelectionOptions atomSelectionOptions;

  @Mixin
  private DynamicsOptions dynamicsOptions;

  @Mixin
  private WriteoutOptions writeOutOptions;

  @Mixin
  private RepExOptions repEx;

  /**
   * --pH or --constantPH Constant pH value for molecular dynamics.
   */
  @Option(names = {"--pH", "--constantPH"}, paramLabel = "7.4", defaultValue = "7.4",
      description = "Constant pH value for molecular dynamics")
  private double pH;

  @Option(names = {"--titrationSteps"}, paramLabel = "10", defaultValue = "10",
      description = "Number of steps done titrating protons on CPU in one cycle")
  private int titrSteps;

  @Option(names = {"--coordinateSteps"}, paramLabel = "100", defaultValue = "100",
      description = "Number of steps done propagating coordinates only on GPU in one cycle")
  private int coordSteps;

  @Option(names = {"--cycles", "--OMMcycles"}, paramLabel = "5", defaultValue = "5",
      description = "Number of times to cycle between titrating protons on CPU and propagating coordinates only on GPU")
  private int cycles;

  @Option(names = {"--titrationReport", "--esvLog"}, paramLabel = "0.25 (psec)", defaultValue = "0.25",
      description = "Interval in psec to report ESV energy and lambdas when cycling between GPU and CPU.")
  private double titrReport;

  @Option(names = {"--pHGaps"}, paramLabel = "1", defaultValue = "1",
      description = "pH gap between replica exchange windows.")
  private double pHGap;

  @Option(names = {"--initDynamics"}, paramLabel = "10000", defaultValue = "10000",
      description = "Number of initialization steps to take before replica exchange windows start.")
  private int initDynamics;

  @Option(names = "--sort", paramLabel = "false", defaultValue = "false",
      description = "Sort archive files by pH")
  private boolean sort;

  @Option(names = "--printRatioData", paramLabel = "false", defaultValue = "false",
      description = "Print out the protonation ratios from throughout the simulation at the end")
  private boolean printRatio;

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = "XYZ or PDB input files.")
  private String filename;

  /**
   * Creation of a public field to try and make the JUnit test work, original code does not declare this as a public field.
   * Originally it is declared in the run method
   */
  private Potential potential;
  public MolecularDynamics molecularDynamics = null;

  /**
   * PhDynamics Constructor.
   */
  public PhDynamics() {
    super();
  }

  /**
   * PhDynamics Constructor.
   * @param binding The Binding to use.
   */
  public PhDynamics(FFXBinding binding) {
    super(binding);
  }

  /**
   * PhDynamics constructor that sets the command line arguments.
   * @param args Command line arguments.
   */
  public PhDynamics(String[] args) {
    super(args);
  }

  public MolecularDynamics getMolecularDynamics() {
    return molecularDynamics;
  }

  public Potential getPotentialObject() {
    return potential;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public PhDynamics run() {

    if (!init()) {
      return this;
    }

    dynamicsOptions.init();

    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }
    // Set the filename.
    String filename = activeAssembly.getFile().getAbsolutePath();
    // Set active atoms.
    atomSelectionOptions.setActiveAtoms(activeAssembly);

    potential = activeAssembly.getPotentialEnergy();

    // Restart File
    File esv = new File(FilenameUtils.removeExtension(filename) + ".esv");
    if (!esv.exists()) {
      esv = null;
    }

    // Initialize and attach extended system first.
    ExtendedSystem esvSystem = new ExtendedSystem(activeAssembly, pH, esv);
    esvSystem.setConstantPh(pH);
    if (potential instanceof ForceFieldEnergy) {
      ((ForceFieldEnergy) potential).attachExtendedSystem(esvSystem);
    }

    int numESVs = esvSystem.getExtendedResidueList().size();
    logger.info(format(" Attached extended system with %d residues.", numESVs));

    double[] x = new double[potential.getNumberOfVariables()];
    potential.getCoordinates(x);
    potential.energy(x, true);
    SystemFilter systemFilter = algorithmFunctions.getFilter();
    if (systemFilter instanceof XYZFilter) {
      XPHFilter xphFilter = new XPHFilter(activeAssembly.getFile(), activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties(), esvSystem);
      xphFilter.readFile();
      logger.info("Reading ESV lambdas from XPH file");
      potential.getCoordinates(x);
      potential.energy(x, true);
    }

    logger.info("\n Running molecular dynamics on " + filename);

    molecularDynamics =
        dynamicsOptions.getDynamics(writeOutOptions, potential, activeAssembly, algorithmListener);

    molecularDynamics.attachExtendedSystem(esvSystem, titrReport);

    // Restart File
    File dyn = new File(FilenameUtils.removeExtension(filename) + ".dyn");
    if (!dyn.exists()) {
      dyn = null;
    }

    if (!(molecularDynamics instanceof MolecularDynamicsOpenMM)) {
      if (repEx.getRepEx()) {
        Comm world = Comm.world();
        int size = world.size();
        if (size == 1) {
          System.exit(0);
        }

        String pHWindows = Arrays.toString(PhReplicaExchange.setEvenSpacePhLadder(pHGap, pH, size));
        pHWindows = pHWindows.replace("[", "")
            .replace("]", "")
            .replace(",", " ");
        CompositeConfiguration properties = activeAssembly.getProperties();
        pHWindows = properties.getString("pH.Windows", pHWindows);
        String[] temp = pHWindows.split(" +");
        if (temp.length != size) {
          logger.severe("pHLadder specified in properties/key file has " +
              "incorrect number of windows given world.size()");
        }
        double[] pHLadder = new double[size];
        for (int i = 0; i < temp.length; i++) {
          pHLadder[i] = Double.parseDouble(temp[i]);
        }

        logger.info("\n Running replica exchange molecular dynamics on " + filename);
        int rank = world.rank();

        File structureFile = new File(filename);
        final String newMolAssemblyFile = structureFile.getParent() + File.separator + rank +
            File.separator + structureFile.getName();
        logger.info(" Set activeAssembly filename: " + newMolAssemblyFile);
        activeAssembly.setFile(new File(newMolAssemblyFile));
        
        try {
          PhReplicaExchange pHReplicaExchange = new PhReplicaExchange(molecularDynamics, structureFile, pH, pHLadder,
              dynamicsOptions.getTemperature(), esvSystem, x, world.size());

          long totalSteps = dynamicsOptions.getNumSteps();
          int nSteps = repEx.getReplicaSteps();
          int exchangeCycles = (int) (totalSteps / nSteps);

          if (exchangeCycles <= 0) {
            exchangeCycles = 1;
          }

          pHReplicaExchange.sample(exchangeCycles, nSteps, dynamicsOptions.getDt(), dynamicsOptions.getReport(),
                  dynamicsOptions.getDt() * titrSteps - 1, initDynamics);

          if (sort) {
            sortMyArc(structureFile, size, pHReplicaExchange.getPhScale()[world.rank()], world.rank());
          }
        } catch (Exception e) {
          logger.severe(" Error during pH replica exchange: " + e.getMessage());
          return this;
        }

      } else {
        // CPU Constant pH Dynamics
        molecularDynamics.dynamic(dynamicsOptions.getSteps(), dynamicsOptions.getDt(),
            dynamicsOptions.getReport(), dynamicsOptions.getWrite(), dynamicsOptions.getTemperature(), true, dyn);
        esvSystem.writeLambdaHistogram(true);
      }
    } else {
      // CPU Constant pH Dynamics alternatives with GPU Dynamics at fixed protonation states.

      // Save a reference to the OpenMM Molecular Dynamics
      MolecularDynamicsOpenMM molecularDynamicsOpenMM = (MolecularDynamicsOpenMM) molecularDynamics;

      // Create an FFX Molecular Dynamics
      molecularDynamics = dynamicsOptions.getDynamics(writeOutOptions, potential,
          activeAssembly, algorithmListener, MDEngine.FFX);
      molecularDynamics.attachExtendedSystem(esvSystem, titrReport);

      if (repEx.getRepEx()) {
        Comm world = Comm.world();
        int size = world.size();

        logger.info("\n Running replica exchange molecular dynamics on " + filename);
        int rank = (size > 1) ? world.rank() : 0;
        logger.info(" Rank: " + rank);

        String pHWindows = Arrays.toString(PhReplicaExchange.setEvenSpacePhLadder(pHGap, pH, size));
        pHWindows = pHWindows.replace("[", "")
            .replace("]", "")
            .replace(",", " ");
        CompositeConfiguration properties = activeAssembly.getProperties();
        pHWindows = properties.getString("pH.Windows", pHWindows);
        String[] temp = pHWindows.split(" +");
        if (temp.length != size) {
          logger.severe("pHLadder specified in properties/key file has " +
              "incorrect number of windows given world.size()");
        }
        double[] pHLadder = new double[size];
        for (int i = 0; i < temp.length; i++) {
          pHLadder[i] = Double.parseDouble(temp[i]);
        }

        File structureFile = new File(filename);
        File rankDirectory = new File(structureFile.getParent() + File.separator + Integer.toString(rank));

        final String newMolAssemblyFile = rankDirectory.getPath() + File.separator + structureFile.getName();
        logger.info(" Set activeAssembly filename: " + newMolAssemblyFile);
        activeAssembly.setFile(new File(newMolAssemblyFile));
        
        try {
          PhReplicaExchange pHReplicaExchange = new PhReplicaExchange(molecularDynamics, structureFile, pH, pHLadder,
              dynamicsOptions.getTemperature(), esvSystem, x, molecularDynamicsOpenMM, potential, world.size());

          pHReplicaExchange.
              sample(cycles, titrSteps, coordSteps, dynamicsOptions.getDt(), dynamicsOptions.getReport(),
                  dynamicsOptions.getDt() * titrSteps, initDynamics);

          if (sort) {
            sortMyArc(structureFile, world.size(), pHReplicaExchange.getPhScale()[world.rank()], world.rank());
          }
        } catch (Exception e) {
          logger.severe("Error during pH replica exchange with OpenMM: " + e.getMessage());
          return this;
        }

      } else {
        for (int i = 0; i < cycles; i++) {
          // Try running on the CPU
          molecularDynamics.setCoordinates(x);
          double forceWriteInterval = titrSteps * 0.001;
          molecularDynamics.dynamic(titrSteps, dynamicsOptions.getDt(), titrReport, forceWriteInterval,
              dynamicsOptions.getTemperature(), true, dyn);
          x = molecularDynamics.getCoordinates();
          esvSystem.writeLambdaHistogram(true);

          // Try running in OpenMM
          potential.energy(x);
          molecularDynamicsOpenMM.setCoordinates(x);
          if (coordSteps != 0) {
            molecularDynamicsOpenMM.dynamic(coordSteps, dynamicsOptions.getDt(), dynamicsOptions.getReport(),
                dynamicsOptions.getWrite(), dynamicsOptions.getTemperature(), true, dyn);
            x = molecularDynamicsOpenMM.getCoordinates();
          }
        }
      }
    }

    if (printRatio) {
      esvSystem.printProtonationRatios();
    }

    return this;
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
      potentials = Collections.singletonList(potential);
    }
    return potentials;
  }

  /**
   * Sort archive files by pH with string parsing
   */
  public static void sortMyArc(File structureFile, int nReplicas, double pH, int myRank) {
    logger.info(" Sorting archive for rank " + myRank);
    String parent = structureFile.getParent();
    String arcName = FilenameUtils.removeExtension(structureFile.getName()) + ".arc";
    BufferedReader[] bufferedReaders = new BufferedReader[nReplicas];
    File output = new File(parent + File.separator + myRank + File.separator + arcName + "_sorted");

    try (BufferedWriter out = new BufferedWriter(new FileWriter(output))) {
      // Get snap length from first directory
      File temp = new File(parent + File.separator + 0 + File.separator + arcName);
      int snapLength = 0;
      int totalLines = 0;
      try (BufferedReader brTemp = new BufferedReader(new FileReader(temp))) {
        String data = brTemp.readLine();
        boolean startSnap = false;
        while (data != null) {
          if (data.contains("pH:")) {
            startSnap = !startSnap;
            if (!startSnap) {
              break;
            }
          }
          data = brTemp.readLine();
          snapLength++;
        }
        totalLines = snapLength;
        while (data != null) {
          totalLines++;
          data = brTemp.readLine();
        }
      }
      int numSnaps = totalLines / snapLength;

      // Build file readers
      for (int i = 0; i < nReplicas; i++) {
        File file = new File(parent + File.separator + i + File.separator + arcName);
        bufferedReaders[i] = new BufferedReader(new FileReader(file));
      }

      try {
        // Read all arc files one snap at a time
        for (int i = 0; i < numSnaps; i++) {
          for (int j = 0; j < nReplicas; j++) {
            // Read up to the first line
            String data = bufferedReaders[j].readLine();
            while (data != null) {
              if (data.contains("pH:")) {
                break;
              }
              data = bufferedReaders[j].readLine();
            }
            // Get pH from line
            String[] tokens = data.split(" +");
            double snapPh = Double.parseDouble(tokens[tokens.length - 3]); // FIXME: pH is not the last index of tokens 'Rank: #' is

            // Add lines to file if correct, otherwise don't
            for (int k = 0; k < snapLength - 1; k++) {
              if (snapPh == pH) {
                out.write(data + "\n");
              }// Readlines doesn't work as expected
              data = bufferedReaders[j].readLine();
            }
            if (snapPh == pH) {
              out.write("\n");
            }
          }
        }
      } catch (IOException e) {
        e.printStackTrace();
      } finally {
        // Cleanup
        for (int i = 0; i < nReplicas; i++) {
          if (bufferedReaders[i] != null) {
            try {
              bufferedReaders[i].close();
            } catch (IOException e) {
              e.printStackTrace();
            }
          }
        }
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
