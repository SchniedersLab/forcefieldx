//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.algorithms.groovy.test

import edu.rit.mp.DoubleBuf
import edu.rit.pj.Comm
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.dynamics.MolecularDynamics
import ffx.algorithms.dynamics.MDEngine
import ffx.algorithms.dynamics.MolecularDynamicsOpenMM
import ffx.numerics.Potential
import ffx.potential.bonded.Residue
import ffx.potential.cli.WriteoutOptions
import ffx.potential.extended.ExtendedSystem
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import java.util.logging.Level

import static java.lang.String.format

/**
 * Use the BAR algorithm to estimate a free energy difference for a CpHMD system.
 * <br>
 * Usage: Umod calculation for model compounds
 * <br>
 * ffxc test.PhBar [options] &lt;filename&gt [file2...];
 */
@Command(description = " Use the BAR algorithm to estimate a free energy difference for a CpHMD system.", name = "test.PhBar")
class PhBar extends AlgorithmsScript {

  @Mixin
  DynamicsOptions dynamicsOptions

  @Mixin
  WriteoutOptions writeOutOptions

  /**
   * --pH or --constantPH Constant pH value for molecular dynamics.
   */
  @Option(names = ['--pH', '--constantPH'], paramLabel = '7.4',
      description = 'Constant pH value for molecular dynamics windows')
  double pH = 7.4

  /**
   * --titrationFix Constant titration state for windows [0-10]
   */
  @Option(names = ['--titrationFix'], paramLabel = '0',
      description = 'Constant titration value for molecular dynamics windows')
  double fixedTitrationState = -1

  /**
   * --tautomerFix Constant tautomer state for windows [0-10]
   */
  @Option(names = ['--tautomerFix'], paramLabel = '0',
      description = 'Constant tautomer value for molecular dynamics windows')
  double fixedTautomerState = -1

  @Option(names = ['--iterations'], paramLabel = '999',
      description = 'Number of times to evaluate neighbor energies')
  int cycles = 999

  @Option(names = ['--coordinateSteps'], paramLabel = '10000',
      description = 'Number of steps done on GPU before each evaluation')
  int coordSteps = 10000

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = "XYZ or PDB input files.")
  private String filename


  private Potential forceFieldEnergy
  MolecularDynamics molecularDynamics = null
  boolean lockTitration = false
  boolean lockTautomer = false
  double energy
  ArrayList<Double> previous = new ArrayList<>()
  ArrayList<Double> current = new ArrayList<>()
  ArrayList<Double> next = new ArrayList<>()


  /**
   * Thermodynamics Constructor.
   */
  PhBar() {
    this(new Binding())
  }

  /**
   * Thermodynamics Constructor.
   * @param binding The Groovy Binding to use.
   */
  PhBar(Binding binding) {
    super(binding)
  }

  PhBar run() {

    if (!init()) {
      return this
    }

    Comm world = Comm.world()
    int nRanks = world.size()
    if (nRanks < 2) {
      logger.severe(" Running BAR with less then the required amount of windows")
    }
    int myRank = (nRanks > 1) ? world.rank() : 0

    // Init titration states and decide what is being locked
    if (fixedTitrationState == -1 && fixedTautomerState == -1) {
      logger.severe(
          " Must select a tautomer or titration to fix windows at. The program will not continue")
      return this
    } else if (fixedTautomerState != -1) {
      fixedTitrationState = (double) myRank / nRanks
      lockTautomer = true
      logger.info(
          " Running BAR across titration states with tautomer state locked at " + fixedTautomerState)
      logger.info(" Titration state for this rank(" + myRank + "): " + fixedTitrationState)
    } else if (fixedTitrationState != -1) {
      fixedTautomerState = (double) myRank / nRanks
      lockTitration = true
      logger.info(" Running BAR across tautomer states with titration state locked at " +
          fixedTitrationState)
      logger.info(" Tautomer state for this rank: " + fixedTautomerState)
    }

    // Illegal state checks
    if (fixedTautomerState > 1 || fixedTitrationState > 1) {
      logger.severe(" ERROR: Cannot assign lambda state to > 1")
      return this
    }

    if (fixedTautomerState < 0 || fixedTitrationState < 0) {
      logger.severe(" ERROR: Cannot assign lambda state to < 1")
      return this
    }

    dynamicsOptions.init()

    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }
    forceFieldEnergy = activeAssembly.getPotentialEnergy()

    // Set the filename.
    String filename = activeAssembly.getFile().getAbsolutePath()

    // Initialize and attach extended system first.
    ExtendedSystem esvSystem = new ExtendedSystem(activeAssembly, pH, null)

    //Setting the systems locked states
    for (Residue res : esvSystem.getExtendedResidueList()) {
      esvSystem.setTitrationLambda(res, fixedTitrationState)
      if (esvSystem.isTautomer(res)) {
        esvSystem.setTautomerLambda(res, fixedTautomerState)
      }
    }

    esvSystem.setConstantPh(pH)
    esvSystem.setFixedTitrationState(true)
    esvSystem.setFixedTautomerState(true)
    forceFieldEnergy.attachExtendedSystem(esvSystem)

    int numESVs = esvSystem.extendedResidueList.size()
    logger.info(format(" Attached extended system with %d residues.", numESVs))

    double[] x = new double[forceFieldEnergy.getNumberOfVariables()]
    forceFieldEnergy.getCoordinates(x)
    forceFieldEnergy.energy(x, true)

    logger.info("\n Running molecular dynamics on " + filename)

    molecularDynamics =
        dynamicsOptions
            .getDynamics(writeOutOptions, forceFieldEnergy, activeAssembly, algorithmListener)

    File structureFile = new File(filename)
    File rankDirectory = new File(
        structureFile.getParent() + File.separator + Integer.toString(myRank))
    if (!rankDirectory.exists()) {
      rankDirectory.mkdir()
    }
    final String newMolAssemblyFile =
        rankDirectory.getPath() + File.separator + structureFile.getName()
    File newDynFile = new File(rankDirectory.getPath() + File.separator +
        FilenameUtils.removeExtension(structureFile.getName()) + ".dyn")
    logger.info(" Set activeAssembly filename: " + newMolAssemblyFile)
    activeAssembly.setFile(new File(newMolAssemblyFile))

    File natNMinusOne = new File(
        rankDirectory.getPath() + File.separator + myRank + "at" + (myRank - 1) + ".log")
    File natN = new File(rankDirectory.getPath() + File.separator + myRank + "at" + myRank + ".log")
    File natNPlusOne = new File(
        rankDirectory.getPath() + File.separator + myRank + "at" + (myRank + 1) + ".log")

    File[] files = new File[] {natNMinusOne, natN, natNPlusOne}
    ArrayList<Double>[] lists = new ArrayList<Double>[] {previous, current, next}

    if (natNMinusOne.exists() && natN.exists() && natNPlusOne.exists()) {
      //MolA should be set to the previous arc file
      logger.warning(
          " Initializing to last run. Since this program deletes these files on completion, runs in which these files already exist are assumed to be restarts.")

      // Read in previous log files into this program
      for (int i = 0; i < 3; i++) {
        try (FileReader fr = new FileReader(files[i])
            BufferedReader br = new BufferedReader(fr)) {
          String data = br.readLine()
          while (data != null) {
            if (myRank != 0 && myRank != world.size() - 1) {
              lists[i].add(Double.parseDouble(data))
              data = br.readLine()
            } else if (myRank == 0) {
              if (i != 0) {
                lists[i].add(Double.parseDouble(data))
                data = br.readLine()
              }
            } else {
              if (i != 2) {
                lists[i].add(Double.parseDouble(data))
                data = br.readLine()
              }
            }
          }
        } catch (IOException e) {
          e.printStackTrace()
        }
      }

      // Check that the program did not end unevenly
      if (myRank != 0 && myRank != world.size() - 1) {
        if (current.size() != next.size() || current.size() != previous.size()) {
          logger.severe(" Restart nAtN.log files are not of the same size")
        } else {
          logger.info(
              " Rank " + myRank + " starting from cycle " + (current.size() + 1) + " of " + cycles)
        }
      } else if (myRank == 0) {
        if (current.size() != next.size()) {
          logger.severe(" Restart nAtN.log files are not of the same size")
        } else {
          logger.info(
              " Rank " + myRank + " starting from cycle " + (current.size() + 1) + " of " + cycles)
        }
      } else {
        if (current.size() != previous.size()) {
          logger.severe(" Restart nAtN.log files are not of the same size")
        } else {
          logger.info(
              " Rank " + myRank + " starting from cycle " + (current.size() + 1) + " of " + cycles)
        }
      }
    }

    if (molecularDynamics instanceof MolecularDynamicsOpenMM) {
      MolecularDynamicsOpenMM molecularDynamicsOpenMM = molecularDynamics

      molecularDynamics =
          dynamicsOptions
              .getDynamics(writeOutOptions, forceFieldEnergy, activeAssembly, algorithmListener,
                  MDEngine.FFX)
      molecularDynamics.attachExtendedSystem(esvSystem, dynamicsOptions.report)

      // Cycle MD runs
      for (int i = current.size(); i < cycles; i++) {
        logger.info(" --------------------------- Cycle " + (i + 1) + " out of " + cycles +
            " ---------------------------")

        molecularDynamics.setCoordinates(x)
        molecularDynamics.dynamic(1, dynamicsOptions.dt, 1, dynamicsOptions.write,
            dynamicsOptions.temperature, true, newDynFile)

        molecularDynamicsOpenMM.setCoordinates(x)
        forceFieldEnergy.energy(x)

        molecularDynamicsOpenMM
            .dynamic(coordSteps, dynamicsOptions.dt, dynamicsOptions.report, dynamicsOptions.write,
                dynamicsOptions.temperature, true, newDynFile)
        x = molecularDynamicsOpenMM.getCoordinates()

        double titrationNeighbor = lockTitration ? 0 : (double) 1 / nRanks
        double tautomerNeighbor = lockTautomer ? 0 : (double) 1 / nRanks

        // Evaluate different energies
        if (myRank != nRanks - 1) {
          for (Residue res : esvSystem.getExtendedResidueList()) {
            esvSystem.setTitrationLambda(res, fixedTitrationState + titrationNeighbor, false)
            //esvSystem.perturbLambdas(lockTautomer, fixedTitrationState + titrationNeighbor) //TODO: get this to work
            if (esvSystem.isTautomer(res)) {
              esvSystem.setTautomerLambda(res, fixedTautomerState + tautomerNeighbor, false)
            }
          }
          forceFieldEnergy.getCoordinates(x)
          energy = forceFieldEnergy.energy(x, false)
          next.add(energy)
        }

        if (myRank != 0) {
          for (Residue res : esvSystem.getExtendedResidueList()) {
            esvSystem.setTitrationLambda(res, fixedTitrationState - titrationNeighbor, false)
            if (esvSystem.isTautomer(res)) {
              esvSystem.setTautomerLambda(res, fixedTautomerState - tautomerNeighbor, false)
            }
          }
          forceFieldEnergy.getCoordinates(x)
          energy = forceFieldEnergy.energy(x, false)
          previous.add(energy)
        }

        for (Residue res : esvSystem.getExtendedResidueList()) {
          esvSystem.setTitrationLambda(res, fixedTitrationState, false)
          if (esvSystem.isTautomer(res)) {
            esvSystem.setTautomerLambda(res, fixedTautomerState, false)
          }
        }
        forceFieldEnergy.getCoordinates(x)
        energy = forceFieldEnergy.energy(x, false)
        current.add(energy)

        logger.info(" Previous: " + previous + "\n Current: " + current + "\n Next: " + next)

        // Write to restart logs
        for (int j = 0; j < 3; j++) {
          try (FileWriter fw = new FileWriter(files[j], true)
              BufferedWriter bw = new BufferedWriter(fw)) {
            if (myRank != 0 && myRank != world.size() - 1) {
              bw.write(lists[j].get(i) as String + "\n")
            } else if (myRank == 0) {
              if (j != 0) {
                bw.write(lists[j].get(i) as String + "\n")
              }
            } else {
              if (j != 2) {
                bw.write(lists[j].get(i) as String + "\n")
              }
            }
          } catch (IOException e) {
            e.printStackTrace()
          }
        }
      }
    } else {
      logger.severe(" MD is not an instance of MDOMM (try adding -Dplatform=OMM --mdE OpenMM)")
    }

    if (current.size() != cycles) {
      logger.severe(" Size of the self energies array is not equal to number of cycles")
      return this
    }

    double[][] parameters = new double[nRanks][current.size() * 3]
    DoubleBuf[] parametersBuf = new DoubleBuf[nRanks]
    int counter = 0
    if (myRank != 0) {
      for (int i = 0; i < previous.size(); i++) {
        parameters[myRank][i] = previous.get(i)
        counter = i
      }
    } else {
      counter = current.size()
    }

    for (int i = 0; i < current.size(); i++) {
      parameters[myRank][counter] = current.get(i)
      counter += 1
    }

    if (myRank != nRanks - 1) {
      for (int i = 0; i < next.size(); i++) {
        parameters[myRank][counter] = next.get(i)
        counter += 1
      }
    }

    for (int i = 0; i < nRanks; i++) {
      parametersBuf[i] = DoubleBuf.buffer(parameters[i])
    }

    DoubleBuf myParametersBuf = parametersBuf[myRank]

    try {
      world.allGather(myParametersBuf, parametersBuf)
    } catch (IOException ex) {
      String message = " CreateBAR allGather failed."
      logger.log(Level.SEVERE, message, ex)
    }

    // File manipulation to create bars
    if (myRank != nRanks - 1) {
      File outputDir = new File(rankDirectory.getParent() + File.separator + "barFiles")
      if (!outputDir.exists()) {
        outputDir.mkdir()
      }
      File output = new File(outputDir.getPath() + File.separator + "energy_" + myRank + ".bar")

      try (FileWriter fw = new FileWriter(output)
          BufferedWriter bw = new BufferedWriter(fw)) {
        bw.write(format("    %d  %f  this.xyz\n", current.size(), dynamicsOptions.temperature))
        for (int i = 1; i <= current.size(); i++) {
          bw.write(format("%5d%17.9f%17.9f\n", i, parameters[myRank][current.size() + i - 1],
              parameters[myRank][current.size() * 2 + i - 2]))
        }

        bw.write(format("    %d  %f  this.xyz\n", current.size(), dynamicsOptions.temperature))
        for (int i = 1; i <= current.size(); i++) {
          bw.write(format("%5d%17.9f%17.9f\n", i, parameters[myRank + 1][i - 1],
              parameters[myRank + 1][current.size() + i - 1]))
        }
      } catch (IOException e) {
        e.printStackTrace()
      }
    }

    return this
  }
}
