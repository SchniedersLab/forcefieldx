//******************************************************************************
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
//******************************************************************************
package ffx.potential.groovy

import com.google.common.collect.MinMaxPriorityQueue
import ffx.crystal.Crystal
import ffx.numerics.Potential
import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The Energy script evaluates the energy of a system.
 * <br>
 * Usage:
 * <br>
 * ffxc Energy &lt;filename&gt;
 */
@Command(description = " Compute the force field potential energy.", name = "ffxc Energy")
class Energy extends PotentialScript {

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  /**
   * -m or --moments print out electrostatic moments.
   */
  @Option(names = ['-m', '--moments'], paramLabel = "false", defaultValue = "false",
      description = 'Print out electrostatic moments.')
  private boolean moments = false

  /**
   * -g or --gradient to print out gradients.
   */
  @Option(names = ['-g', '--gradient'], paramLabel = "false", defaultValue = "false",
      description = 'Compute the atomic gradient as well as energy.')
  private boolean gradient = false

  /**
   * * --fl or --findLowest Return the n lowest energy structures from an ARC or PDB file.
   */
  @Option(names = ['--fl', '--findLowest'], paramLabel = "0", defaultValue = "0",
      description = 'Return the n lowest energies from an ARC/PDB file.')
  private int fl = 0

  /**
   * -v or --verbose enables printing out all energy components for multi-snapshot files (
   * the first snapshot is always printed verbosely).
   */
  @Option(names = ['-v', '--verbose'], paramLabel = "false", defaultValue = "false",
      description = "Print out all energy components for each snapshot.")
  private boolean verbose = false

  /**
   * --dc or --densityCutoff Collect structures above a specified density.
   */
  @Option(names = ['--dc', '--densityCutoff'], paramLabel = "0.0", defaultValue = "0.0",
          description = "Create ARC file of structures above a specified density.")
  private double dCutoff = 0.0

  /**
   * --ec or --energyCutOff Collect structures below a specified energy range from the minimum energy.
   */
  @Option(names = ['--ec', '--energyCutoff'], paramLabel = "0.0", defaultValue = "0.0",
          description = "Create ARC file of structures within a specified energy of the lowest energy structure.")
  private double eCutoff = 0.0

  /**
   * The final argument is a PDB or XYZ coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in PDB or XYZ format.')
  private String filename = null

  public double energy = 0.0
  public ForceFieldEnergy forceFieldEnergy = null
  private AssemblyState assemblyState = null

  /**
   * Energy constructor.
   */
  Energy() {
    this(new Binding())
  }

  /**
   * Energy constructor.
   * @param binding The Groovy Binding to use.
   */
  Energy(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  Energy run() {
    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    logger.info("\n Running Energy on " + filename)

    // Apply atom selections
    atomSelectionOptions.setActiveAtoms(activeAssembly)

    forceFieldEnergy = activeAssembly.getPotentialEnergy()
    int nVars = forceFieldEnergy.getNumberOfVariables()
    double[] x = new double[nVars]
    forceFieldEnergy.getCoordinates(x)

    if (gradient) {
      double[] g = new double[nVars]
      int nAts = (int) (nVars / 3)
      energy = forceFieldEnergy.energyAndGradient(x, g, true)
      logger.info(format("    Atom       X, Y and Z Gradient Components (kcal/mol/A)"))
      for (int i = 0; i < nAts; i++) {
        int i3 = 3 * i
        logger.info(format(" %7d %16.8f %16.8f %16.8f", i + 1, g[i3], g[i3 + 1], g[i3 + 2]))
      }
    } else {
      energy = forceFieldEnergy.energy(x, true)
    }

    if (moments) {
      Atom[] activeAtoms = activeAssembly.getActiveAtomArray()
      forceFieldEnergy.getPmeNode().computeMoments(activeAtoms, false)
    }

    SystemFilter systemFilter = potentialFunctions.getFilter()

    if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {

      int numSnaps = fl
      double lowestEnergy = energy
      assemblyState = new AssemblyState(activeAssembly)
      int index = 1

      // Making the MinMax priority queue that will expel the largest entry when it reaches its maximum size N/
      MinMaxPriorityQueue<StateContainer> lowestEnergyQueue = null
      if (fl > 0) {
        lowestEnergyQueue = MinMaxPriorityQueue.maximumSize(numSnaps).create()
        lowestEnergyQueue.add(new StateContainer(assemblyState, lowestEnergy))
      }

      int numModels = systemFilter.countNumModels()
      //Store densities in ordered encountered (used in density cutoff implementation).
      double[] densities = new double[numModels]
      //Store energies in ordered encountered (used in energy cutoff implementation).
      double[] energies = new double[numModels]
      energies[0] = forceFieldEnergy.energy(x, false)

      while (systemFilter.readNext()) {
        index++
        Crystal crystal = activeAssembly.getCrystal()
        densities[index - 1] = crystal.getDensity(activeAssembly.getMass())
        forceFieldEnergy.setCrystal(crystal)
        forceFieldEnergy.getCoordinates(x)
        if (verbose) {
          logger.info(format(" Snapshot %4d", index))
          if(!crystal.aperiodic()){
            logger.info(format("\n Density:                                %6.3f (g/cc)",
                    crystal.getDensity(activeAssembly.getMass())))
          }
          energy = forceFieldEnergy.energy(x, true)
        } else {
          energy = forceFieldEnergy.energy(x, false)
          logger.info(format(" Snapshot %4d: %16.8f (kcal/mol)", index, energy))
        }

        energies[index - 1] = energy

        //Update lowest encountered energy
        if (energy < lowestEnergy) {
          lowestEnergy = energy
        }

        if (fl > 0) {
          lowestEnergyQueue.add(new StateContainer(new AssemblyState(activeAssembly), energy))
        }

        if (moments) {
          Atom[] activeAtoms = activeAssembly.getActiveAtomArray()
          forceFieldEnergy.getPmeNode().computeMoments(activeAtoms, false)
        }

      }

      // If cutoffs have been selected create an ARC or PDB to store structures that satisfy cutoff
      if((eCutoff>0.0 || dCutoff>0.0) && numModels > 1) {
        systemFilter.readNext(true)
        activeAssembly = systemFilter.getActiveMolecularSystem()
        String modelFilename = filenames.get(0)
        String fileName = FilenameUtils.getName(modelFilename)
        String ext = FilenameUtils.getExtension(fileName)
        fileName = FilenameUtils.removeExtension(fileName)
        File saveFile
        SystemFilter writeFilter
        if (ext.toUpperCase().contains("XYZ")) {
          saveFile = new File(fileName + ".xyz")
          writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
                  activeAssembly.getProperties())
          potentialFunctions.saveAsXYZ(activeAssembly, saveFile)
        } else if (ext.toUpperCase().contains("ARC")) {
          saveFile = new File(fileName + ".arc")
          saveFile = potentialFunctions.versionFile(saveFile)
          writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
                  activeAssembly.getProperties())
          logger.info("SaveFile: " + saveFile.getAbsolutePath())
          saveFile.createNewFile()
        } else {
          saveFile = new File(fileName + ".pdb")
          saveFile = potentialFunctions.versionFile(saveFile)
          writeFilter = new PDBFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
                  activeAssembly.getProperties())
          if (numModels > 1) {
            writeFilter.setModelNumbering(0)
          }
          saveFile.createNewFile()
        }

        // Determine if each structure meets the cutoff condition
        if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {
          int structNum = 0
          if (eCutoff > 0.0) {
            logger.info(format("Lowest Energy of: %16.4f\n Saving structures with energy below: %16.4f", lowestEnergy,
                    lowestEnergy + eCutoff))
          }

          do {
            if (dCutoff > 0.0 && eCutoff > 0.0) {
              if(energies[structNum] < lowestEnergy + eCutoff && densities[structNum] > dCutoff) {
                if (systemFilter instanceof PDBFilter) {
                  PDBFilter pdbFilter = (PDBFilter) systemFilter
                  pdbFilter.writeFile(saveFile, true, false, false)
                  saveFile.append("ENDMDL\n")
                } else if (systemFilter instanceof XYZFilter) {
                  writeFilter.writeFile(saveFile, true)
                }
              }
            } else if (dCutoff > 0.0) {
              if(densities[structNum] > dCutoff) {
                if (systemFilter instanceof PDBFilter) {
                  PDBFilter pdbFilter = (PDBFilter) systemFilter
                  pdbFilter.writeFile(saveFile, true, false, false)
                  saveFile.append("ENDMDL\n")
                } else if (systemFilter instanceof XYZFilter) {
                  writeFilter.writeFile(saveFile, true)
                }
              }
            } else if (eCutoff > 0.0) {
              if(energies[structNum] < lowestEnergy + eCutoff) {
                if (systemFilter instanceof PDBFilter) {
                  PDBFilter pdbFilter = (PDBFilter) systemFilter
                  pdbFilter.writeFile(saveFile, true, false, false)
                  saveFile.append("ENDMDL\n")
                } else if (systemFilter instanceof XYZFilter) {
                  writeFilter.writeFile(saveFile, true)
                }
              }
            }
            structNum++
          } while (systemFilter.readNext())
          if (systemFilter instanceof PDBFilter) {
            saveFile.append("END\n")
          }
        }
      }

      if (fl > 0) {
        if (numSnaps > index) {
          logger.warning(format(
              " Requested %d snapshots, but file %s has only %d snapshots. All %d energies will be reported",
              numSnaps, filename, index, index))
          numSnaps = index
        }

        // Configure the base directory if it has not been set.
        File saveDir = baseDir
        String modelFilename = activeAssembly.getFile().getAbsolutePath()
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
          saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }
        String dirName = saveDir.toString() + File.separator
        String fileName = FilenameUtils.getName(modelFilename)

        for (int i = 0; i < numSnaps - 1; i++) {
          StateContainer savedState = lowestEnergyQueue.removeLast()
          AssemblyState finalAssembly = savedState.getState()
          finalAssembly.revertState()
          double finalEnergy = savedState.getEnergy()
          logger.info(format(" The potential energy found is %16.8f (kcal/mol)", finalEnergy))
          File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
          MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly()
          potentialFunctions.saveAsPDB(molecularAssembly, saveFile)
        }

        StateContainer savedState = lowestEnergyQueue.removeLast()
        AssemblyState lowestAssembly = savedState.getState()
        lowestEnergy = savedState.getEnergy()

        assemblyState.revertState()
        logger.info(format(" The lowest potential energy found is %16.8f (kcal/mol)", lowestEnergy))

        // Prints our final energy (which will be the lowest energy
        File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
        MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly()
        potentialFunctions.saveAsPDB(molecularAssembly, saveFile)
      }
    }

    return this
  }

  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (forceFieldEnergy == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList(forceFieldEnergy)
    }
    return potentials
  }

  private class StateContainer implements Comparable<StateContainer> {

    private final AssemblyState state
    private final double e

    StateContainer(AssemblyState state, double e) {
      this.state = state
      this.e = e
    }

    AssemblyState getState() {
      return state
    }

    double getEnergy() {
      return e
    }

    @Override
    int compareTo(StateContainer o) {
      return Double.compare(e, o.getEnergy())
    }
  }

}

