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
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XPHFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import java.util.logging.Level
import java.util.logging.LogRecord

import static java.lang.String.format

/**
 * Use the BAR algorithm to estimate a free energy difference for a CpHMD system.
 * <br>
 * Usage: Umod calculation for model compounds
 * <br>
 * ffxc test.RaoBlackwellEstimator [options] &lt;filename&gt [file2...];
 */
@Command(description = " Use the Rao-Blackwell estimator to get a free energy difference for residues in a CpHMD system.", name = "test.RaoBlackwellEstimator")
class RaoBlackwellEstimator extends AlgorithmsScript {

  @Mixin
  DynamicsOptions dynamicsOptions

  @Mixin
  WriteoutOptions writeOutOptions

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = "PDB input file in the same directory as the ARC file.")
  private String filename

  private Potential forceFieldEnergy
  MolecularDynamics molecularDynamics = null
  double energy
  ArrayList<Double> previous = new ArrayList<>()
  ArrayList<Double> current = new ArrayList<>()
  ArrayList<Double> next = new ArrayList<>()

  /**
   * Thermodynamics Constructor.
   */
  RaoBlackwellEstimator() {
    this(new Binding())
  }

  /**
   * Thermodynamics Constructor.
   * @param binding The Groovy Binding to use.
   */
  RaoBlackwellEstimator(Binding binding) {
    super(binding)
  }

  RaoBlackwellEstimator run() {

    if (!init()) {
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
    ExtendedSystem esvSystem = new ExtendedSystem(activeAssembly, 7.0, null)

    // Set up the XPHFilter.
    XPHFilter xphFilter = new XPHFilter(
            activeAssembly.getArchiveFile(),
            activeAssembly,
            activeAssembly.getForceField(),
            activeAssembly.getProperties(),
            esvSystem)

    esvSystem.setFixedTitrationState(true)
    esvSystem.setFixedTautomerState(true)
    forceFieldEnergy.attachExtendedSystem(esvSystem)

    int numESVs = esvSystem.extendedResidueList.size()
    logger.info(format(" Attached extended system with %d residues.", numESVs))

    double[] x = new double[forceFieldEnergy.getNumberOfVariables()]
    forceFieldEnergy.getCoordinates(x)
    forceFieldEnergy.energy(x, true)

    File[] files = new File[] {natNMinusOne, natN, natNPlusOne}
    ArrayList<Double>[] lists = new ArrayList<Double>[] {previous, current, next}

    if (natNMinusOne.exists() && natN.exists() && natNPlusOne.exists()) {
      //MolA should be set to the previous arc file
      logger.warning(
          " Initializing to last run.")

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

        if(i == current.size()) {
          logger.info(" Previous: " + previous + "\n Current: " + current + "\n Next: " + next)
        }

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
    } else if (current.size() == 0 && !createBar){
      logger.severe(" MD is not an instance of MDOMM (try adding -Dplatform=OMM --mdE OpenMM)")
    }

    double[][] parameters = new double[nRanks][cycles * 3]
    DoubleBuf[] parametersBuf = new DoubleBuf[nRanks]
    if (myRank != 0) {
      if(previous.size() != cycles) {
        logger.severe(" Size of the previous energies array is not equal to number of cycles")
        return this
      }
      for (int i = 0; i < previous.size(); i++) {
        parameters[myRank][i] = previous.get(i)
      }
    }

    if(current.size() != cycles) {
      logger.severe(" Size of the current energies array is not equal to number of cycles")
      return this
    }
    for (int i = 0; i < current.size(); i++) {
      parameters[myRank][current.size() + i] = current.get(i)
    }

    if (myRank != nRanks - 1) {
      if(next.size() != cycles) {
        logger.severe(" Size of the next energies array is not equal to number of cycles")
        return this
      }
      for (int i = 0; i < next.size(); i++) {
        parameters[myRank][2 * current.size() + i] = next.get(i)
      }
    }

    logger.info("Parameters: " + parameters[myRank] as String)

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
        for (int i = 0; i < current.size(); i++) {
          bw.write(format("%5d%17.9f%17.9f\n", i+1, parameters[myRank][current.size() + i],
              parameters[myRank][current.size() * 2 + i]))
        }

        bw.write(format("    %d  %f  this.xyz\n", current.size(), dynamicsOptions.temperature))
        for (int i = 0; i <= current.size(); i++) {
          bw.write(format("%5d%17.9f%17.9f\n", i+1, parameters[myRank + 1][i],
              parameters[myRank + 1][current.size() + i]))
        }
      } catch (IOException e) {
        e.printStackTrace()
      }
    }

    return this
  }
}
