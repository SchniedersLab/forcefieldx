//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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

import edu.rit.pj.Comm
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.cli.RepExOptions
import ffx.algorithms.dynamics.MolecularDynamics
import ffx.algorithms.dynamics.MolecularDynamicsOpenMM
import ffx.algorithms.dynamics.PhReplicaExchange
import ffx.numerics.Potential
import ffx.potential.cli.WriteoutOptions
import ffx.potential.extended.ExtendedSystem
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The PhDynamics script implements constant pH molecular dynamics .
 * <br>
 * Usage:
 * <br>
 * ffxc PhDynamics [options] &lt;filename&gt; [file2...]
 */
@Command(description = " Run constant pH dynamics on a system.", name = "ffxc PHDynamics")
class PhDynamics extends AlgorithmsScript {

  @Mixin
  DynamicsOptions dynamicsOptions

  @Mixin
  WriteoutOptions writeOutOptions

  @Mixin
  RepExOptions repEx

  /**
   * --pH or --constantPH Constant pH value for molecular dynamics.
   */
  @Option(names = ['--pH', '--constantPH'], paramLabel = '7.4',
      description = 'Constant pH value for molecular dynamics')
  double pH = 7.4

  @Option(names = ['--titrationSteps'], paramLabel = '10',
          description = 'Number of steps done titrating protons on CPU in one cycle')
  int titrSteps  = 10

  @Option(names = ['--coordinateSteps'], paramLabel = '100',
          description = 'Number of steps done propagating coordinates only on GPU in one cycle')
  int coordSteps  = 10

  @Option(names = ['--cycles', '--OMMcycles'], paramLabel = '5',
          description = 'Number of times to cycle between titrating protons on CPU and propagating coordinates only on GPU')
  int cycles  = 5

  @Option(names = ['--titrationReport', '--esvLog'], paramLabel = '0.001 (psec)',
          description = 'Interval in psec to report ESV energy and lambdas when cycling between GPU and CPU.')
  double titrReport  = 0.001

  @Option(names = ['--pHMin'], paramLabel = '2',
          description = 'Minimum pH for replica exchange windows.')
  double pHMin = 2

  @Option(names = ['--pHMax'], paramLabel = '12',
          description = 'Maximum pH for replica exchange windows.')
  double pHMax = 12



  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = "XYZ or PDB input files.")
  private String filename

  /**
   * Creation of a public field to try and make the JUnit test work, original code does not declare this as a public field.
   * Originally it is declared in the run method
   */
  private Potential potential
  public MolecularDynamics molecularDynamics = null

  MolecularDynamics getMolecularDynamics() {
    return molecularDynamics
  }

  Potential getPotentialObject() {
    return potential
  }

  /**
   * Dynamics Constructor.
   */
  PhDynamics() {
    this(new Binding())
  }

  /**
   * Dynamics Constructor.
   * @param binding The Groovy Binding to use.
   */
  PhDynamics(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  PhDynamics run() {

    if (!init()) {
      return this
    }

    dynamicsOptions.init()

    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }
    potential = activeAssembly.getPotentialEnergy()
    // Set the filename.
    String filename = activeAssembly.getFile().getAbsolutePath()

    // Restart File
    File esv = new File(FilenameUtils.removeExtension(filename) + ".esv")
    if (!esv.exists()) {
      esv = null
    }
    // Initialize and attach extended system first.
    ExtendedSystem esvSystem = new ExtendedSystem(activeAssembly, esv)
    esvSystem.setConstantPh(pH)
    potential.attachExtendedSystem(esvSystem)

    int numESVs = esvSystem.extendedResidueList.size()
    logger.info(format(" Attached extended system with %d residues.", numESVs))

    double[] x = new double[potential.getNumberOfVariables()]
    potential.getCoordinates(x)
    potential.energy(x, true)

    logger.info("\n Running molecular dynamics on " + filename)

    molecularDynamics =
        dynamicsOptions.getDynamics(writeOutOptions, potential, activeAssembly, algorithmListener)

    molecularDynamics.attachExtendedSystem(esvSystem, dynamicsOptions.report)

    // Restart File
    File dyn = new File(FilenameUtils.removeExtension(filename) + ".dyn")
    if (!dyn.exists()) {
      dyn = null
    }


    //TODO: **************************************************************************************
    if(repEx.repEx){
      Comm world = Comm.world()
      int size = world.size()

      logger.info("\n Running replica exchange molecular dynamics on " + filename)
      int rank = (size > 1) ? world.rank() : 0
      logger.info("Rank: " + rank.toString())

      File structureFile = new File(filename)
      File rankDirectory = new File(structureFile.getParent() + File.separator + Integer.toString(rank))
      if (!rankDirectory.exists()) {
        rankDirectory.mkdir()
      }

      final String newMolAssemblyFile = rankDirectory.getPath() + File.separator + structureFile.getName()
      logger.info("Set activeAssembly filename: " + newMolAssemblyFile)
      activeAssembly.setFile(new File(newMolAssemblyFile))
      PhReplicaExchange pHReplicaExchange = new PhReplicaExchange(molecularDynamics, pH, pHMin, pHMax, dynamicsOptions.temperature, esvSystem)

      long totalSteps = dynamicsOptions.numSteps
      int nSteps = repEx.replicaSteps
      int exchangeCycles = (int) (totalSteps / nSteps)

      if (exchangeCycles <= 0) {
        exchangeCycles = 1
      }

      pHReplicaExchange.
              sample(exchangeCycles, nSteps, dynamicsOptions.dt, dynamicsOptions.report, dynamicsOptions.write)


      //TODO: Output xyz files into proper rank directories
      String xyzName
      xyzName = structureFile.getParent() + File.separator + FilenameUtils.removeExtension(structureFile.getName()) + ".xyz"

      String outputName = rankDirectory.getPath() + File.separator + "rankOutput.log"
      File output = new File(outputName)
      try(FileReader r = new FileReader(structureFile.getParent() + File.separator + "repEx.log")
          BufferedReader br = new BufferedReader(r)

          FileWriter wr = new FileWriter(output)
          BufferedWriter bwr = new BufferedWriter(wr)) {
        bwr.write("")
        String data = br.readLine()
        while(data != null) {
          if (data.contains("[" + world.rank() + "]")) {
            wr.write(data + "\n")
          }
          data = br.readLine()
        }
      } catch(IOException e){
        e.printStackTrace()
      }
    } else if (!(molecularDynamics instanceof MolecularDynamicsOpenMM)) {
      // CPU Constant pH Dynamics
      molecularDynamics.dynamic(dynamicsOptions.steps, dynamicsOptions.dt,
          dynamicsOptions.report, dynamicsOptions.write, dynamicsOptions.temperature, true, dyn)
      esvSystem.writeLambdaHistogram()
    } else {
      // CPU Constant pH Dynamics alternatives with GPU Dynamics at fixed protonation states.

      // Save a reference to the OpenMM Molecular Dynamics
      MolecularDynamicsOpenMM molecularDynamicsOpenMM = molecularDynamics

      // Create an FFX Molecular Dynamics
      molecularDynamics =
          dynamicsOptions.getDynamics(writeOutOptions, potential, activeAssembly, algorithmListener,
              MolecularDynamics.DynamicsEngine.FFX)
      molecularDynamics.attachExtendedSystem(esvSystem, titrReport)

      for (int i = 0; i < cycles; i++) {
        // Try running on the CPU
        molecularDynamics.setCoordinates(x)
        molecularDynamics.dynamic(titrSteps, dynamicsOptions.dt, titrReport, dynamicsOptions.write,
                dynamicsOptions.temperature, true, dyn)
        x = molecularDynamics.getCoordinates()
        esvSystem.writeLambdaHistogram()

        // Try running in OpenMM
        potential.energy(x)
	    molecularDynamicsOpenMM.setCoordinates(x)
        if(coordSteps != 0){
          molecularDynamicsOpenMM.dynamic(coordSteps, dynamicsOptions.dt, dynamicsOptions.report, dynamicsOptions.write,
                  dynamicsOptions.temperature, true, dyn)
          x = molecularDynamicsOpenMM.getCoordinates()
        }
      }
    }

    return this
  }

  /**
   * {@inheritDoc}
   */
  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (potential == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList(potential)
    }
    return potentials
  }

}
