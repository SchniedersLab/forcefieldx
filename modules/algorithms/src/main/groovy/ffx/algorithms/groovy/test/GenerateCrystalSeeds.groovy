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

import com.google.common.collect.MinMaxPriorityQueue
import ffx.algorithms.AlgorithmListener
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.optimize.Minimize
import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.RestraintBond
import ffx.potential.nonbonded.pme.Polarization
import ffx.potential.parameters.BondType
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Molecule
import ffx.potential.bonded.Bond

import static ffx.potential.utils.Superpose.applyRotation
import static org.apache.commons.math3.util.FastMath.cos
import static org.apache.commons.math3.util.FastMath.sin

/**
 * GenerateCrystalSeeds is a Groovy script that generates a set of molecular orientations in vacuum and
 * calculates the energy of each conformation.
 * <br>
 * Usage:
 * <br>
 * ffxc test.GenerateCrystalSeeds [options] &lt;filename&gt;
 */
@Command(description = " Calculates interaction energies of different molecular orientations and saves low energy orientations.",
        name = "test.GenerateCrystalSeeds")
class GenerateCrystalSeeds extends AlgorithmsScript {

    /**
     * --eps
     */
    @Option(names = ['--eps'], paramLabel = '.1',
            description = 'Cutoff for minimization.')
    double eps = 0.1

    /**
     * --maxIter
     */
    @Option(names = ['--maxIter'], paramLabel = '100',
            description = 'Max iterations for minimization.')
    int maxIter =  1000

    /**
     * --hBondDist
     */
    @Option(names = ['--hBondDist', '--hbd'], paramLabel = '2.0',
            description = 'Initial h-bond distance.')
    double hBondDist = 2.0

    /**
     * --saveEnergyCutoff
     */
    @Option(names = ['--saveEnergyCutoff', "--sec"], paramLabel = '1.5',
            description = 'Cutoff conformations that have an energy (kcal/mol) less then this cutoff.')
    double saveEnergyCutoff = 1.5

    /**
     * --flatBottomRadius
     */
    @Option(names = ['--flatBottomRadius', "--fbr"], paramLabel = '0.5',
            description = 'Radius of flat bottom potential.')
    double flatBottomRadius = 0.5

    /**
     * --saveNumStates
     */
    @Option(names = ['--saveNumStates', "--sns"], paramLabel = '-1',
            description = 'Save this many of the lowest energy states. This is not the default mode.')
    int saveNumStates = -1

    /**
     * --saveAll
     */
    @Option(names = ['--saveAll', "--sa"], paramLabel = 'true',
            description = 'Default mode.')
    boolean saveAll = true

    /**
     * --excludeH
     */
    @Option(names = ['--excludeH', "--eh"], paramLabel = 'false', defaultValue = 'false',
            description = 'Only include H bonded to electronegative atoms in conformations.')
    private boolean excludeH = false

    /**
     * One filename.
     */
    @Parameters(arity = "1..*", paramLabel = "files",
            description = "XYZ input file.")
    private String filename

    /**
    * Constructor.
    */
  GenerateCrystalSeeds() {
    this(new Binding())
  }

  /**
    * Constructor.
    * @param binding The Groovy Binding to use.
    */
    GenerateCrystalSeeds(Binding binding) {
    super(binding)
  }

  /**
    * {@inheritDoc}
    */
  @Override
  GenerateCrystalSeeds run() {
      // Set system properties as soon as script starts
      // Init the context and bind variables.
      if (!init()) {
          return this
      }


      // Load the MolecularAssembly of the input file.
      activeAssembly = getActiveAssembly(filename)
      if (activeAssembly == null) {
          logger.info(helpString())
          return this
      }

      // Set the filename.
      filename = activeAssembly.getFile().getAbsolutePath()
      activeAssembly.update()

      // Get number of molecules in file
      int numMolecules = activeAssembly.getMoleculeArray().length
      logger.info("\n Found " + numMolecules + " molecules.")
      if (numMolecules != 2){
          logger.severe(" XYZ file must contain 2 molecules")
      }

      // Get full atom lists of both molecules to act (rotation/translation) on
      Molecule[] molecules = activeAssembly.getMoleculeArray()
      Atom[] moleculeOneAtoms = molecules[0].getAtomList()
      Atom[] moleculeTwoAtoms = molecules[1].getAtomList()

      // Add atoms to moleculeOneHeavyAtoms and moleculeTwoHeavyAtoms from filteredAtoms
      ArrayList<Atom> moleculeOneFiltered = new ArrayList<>()
      ArrayList<Atom> moleculeTwoFiltered = new ArrayList<>()
      Atom[] atoms = activeAssembly.getAtomList()
      for(Atom a: atoms){
          if(a.getAtomType().atomicNumber == 7|| a.getAtomType().atomicNumber == 8
                  || a.getAtomType().atomicNumber == 9 || a.getAtomType().atomicNumber == 15
          || a.getAtomType().atomicNumber == 16 || a.getAtomType().atomicNumber == 17){ // NOFPSCl
              if(a.getMoleculeNumber() == activeAssembly.getMoleculeNumbers()[0]) {
                  moleculeOneFiltered.add(a)
              } else{
                  moleculeTwoFiltered.add(a)
              }

              // Searching for bonded h's --> program runs fast enough that we can afford to do all
              if(excludeH) {
                  for (Bond b : a.getBonds()) {
                      int num = b.get1_2(a).getAtomType().atomicNumber
                      if (num == 1) {
                          if (a.getMoleculeNumber() == activeAssembly.getMoleculeNumbers()[0]) {
                              moleculeOneFiltered.add(a)
                          } else {
                              moleculeTwoFiltered.add(a)
                          }
                      }
                  }
              }
          }
          else if(a.getAtomType().atomicNumber == 1 && !excludeH){ // Add all H's --> the increased sampling is worth it?
              if(a.getMoleculeNumber() == activeAssembly.getMoleculeNumbers()[0]) {
                  moleculeOneFiltered.add(a)
              } else{
                  moleculeTwoFiltered.add(a)
              }
          }
      }

      ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
      double[] x = new double[forceFieldEnergy.getNumberOfVariables()]
      forceFieldEnergy.getCoordinates(x)
      forceFieldEnergy.energy(x, false)

      // Minimize Monomers, this step is very important since it affects all energies
      // Monomer one energy
      for(Atom a: moleculeTwoAtoms){ a.setUse(false) }
      Minimize monomerMinEngine
      logger.info("\n --------- Minimize Monomer 1 --------- ")
      forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE)
      monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
      monomerMinEngine.minimize(eps, maxIter)
      forceFieldEnergy.getPmeNode().setPolarization(Polarization.DIRECT)
      monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
      monomerMinEngine.minimize(eps, maxIter)
      forceFieldEnergy.getPmeNode().setPolarization(Polarization.MUTUAL)
      monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
      monomerMinEngine.minimize(eps, maxIter)
      for(Atom a: moleculeTwoAtoms){ a.setUse(true) }

      // Monomer two energy
      for(Atom a: moleculeOneAtoms){ a.setUse(false) }
      logger.info("\n --------- Minimize Monomer 2 --------- ")
      forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE)
      monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
      monomerMinEngine.minimize(eps, maxIter)
      forceFieldEnergy.getPmeNode().setPolarization(Polarization.DIRECT)
      monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
      monomerMinEngine.minimize(eps, maxIter)
      forceFieldEnergy.getPmeNode().setPolarization(Polarization.MUTUAL)
      monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
      monomerMinEngine.minimize(eps, maxIter)
      for(Atom a: moleculeOneAtoms){ a.setUse(true) }

      // Minimize Dimer Structure so energy runs on it
      logger.info("\n --------- Minimize Dimer --------- ")
      forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE)
      monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
      monomerMinEngine.minimize(1.0, maxIter)
      forceFieldEnergy.getPmeNode().setPolarization(Polarization.DIRECT)
      monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
      monomerMinEngine.minimize(1.0, maxIter)
      forceFieldEnergy.getPmeNode().setPolarization(Polarization.MUTUAL)
      monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
      monomerMinEngine.minimize(1.0, maxIter)

      logger.info("\n --------- Monomer 1 Energy Breakdown --------- ")
      for(Atom a: moleculeTwoAtoms){ a.setUse(false) }
      forceFieldEnergy.getCoordinates(x)
      double monomerEnergy = forceFieldEnergy.energy(x, true)
      for(Atom a: moleculeTwoAtoms){ a.setUse(true) }


      logger.info("\n --------- Monomer 2 Energy Breakdown --------- ")
      for(Atom a: moleculeOneAtoms){ a.setUse(false) }
      forceFieldEnergy.getCoordinates(x)
      double monomerEnergy2 = forceFieldEnergy.energy(x, true)

      // Get energies logged for init structure
      logger.info("\n --------- Complex Structure Energy Breakdown --------- ")
      forceFieldEnergy.getCoordinates(x)
      double dimerEnergy = forceFieldEnergy.energy(x, true)

      // Log potentials
      logger.info("\n Monomer energy 1: " + monomerEnergy + " kcal/mol")
      monomerEnergy += monomerEnergy2
      logger.info(" Monomer energy 2: " + monomerEnergy2 + " kcal/mol")
      logger.info(" Total of monomer energies: " + monomerEnergy + " kcal/mol")
      logger.info(" Complex energy: " + dimerEnergy + " kcal/mol")

      // Making the MinMax priority queue that will expel the largest entry when it reaches its maximum size
      MinMaxPriorityQueue<StateContainer> lowestEnergyQueue = MinMaxPriorityQueue.
              maximumSize(moleculeOneFiltered.size() * moleculeTwoFiltered.size()).create()
      if(saveAll){
          logger.info("\n Saving all conformations.")
      } else if (saveNumStates != -1 && saveNumStates > 0) {
          lowestEnergyQueue = MinMaxPriorityQueue.maximumSize(saveNumStates).create()
          logger.info("\n Saving the " + saveNumStates + " lowest energy conformations.")
      } else {
          saveEnergyCutoff = Math.abs(saveEnergyCutoff)
          logger.info("\n Saving conformations within " + saveEnergyCutoff + " kcal/mol of lowest energy conformation")
      }

      //Energy Container to store minimized energies
      ArrayList<Double> energies = new ArrayList<>()
      double[] zAxis = new double[]{0,0,1}

      // Loop through interactions between the two molecules --> Not necessarily symmetric but close?
      int loopCounter = 1
      for(Atom a: moleculeOneFiltered){
          // Get center of mass of moleculeOneAtoms
          double[] moleculeOneCOM = getCOM(molecules[0])
          for(int i = 0; i < moleculeOneAtoms.length; i++){
              moleculeOneAtoms[i].move(-moleculeOneCOM[0], -moleculeOneCOM[1], -moleculeOneCOM[2])
          }
          // Get coordinates of a
          double[] aCoords = a.getXYZ().copy().get()
          // Get rotation matrix to align dipole moment with z-axis
          double[][] rotation = getRotationBetween(aCoords, zAxis)
          // Get a reference to the moleculeOneAtoms
          double[] moleculeOneAtomPositions = new double[moleculeOneAtoms.length * 3]
          for(int i = 0; i < moleculeOneAtoms.length; i++){
              moleculeOneAtomPositions[i*3] = moleculeOneAtoms[i].getX()
              moleculeOneAtomPositions[i*3 + 1] = moleculeOneAtoms[i].getY()
              moleculeOneAtomPositions[i*3 + 2] = moleculeOneAtoms[i].getZ()
          }
          applyRotation(moleculeOneAtomPositions, rotation)
          // Move atoms into rotated positions
          for(int i = 0; i < moleculeOneAtoms.length; i++){
              moleculeOneAtoms[i].setXYZ(moleculeOneAtomPositions[i*3],
                      moleculeOneAtomPositions[i*3 + 1],
                      moleculeOneAtomPositions[i*3 + 2])
          }
          // Move moleculeTwoHeavyAtoms atom to origin and rotate it to align with z-axis
          zAxis[2] = -1 // Align with opposite side of axis --> Reset after loop
          for(Atom b: moleculeTwoFiltered){
              logger.info(" Trial " + loopCounter + " out of " +
                      (moleculeTwoFiltered.size() * moleculeOneFiltered.size()))
              loopCounter++
              // Get center of mass of moleculeTwoAtoms and set to zero
              double[] moleculeTwoCOM = getCOM(molecules[1])
              for(int i = 0; i < moleculeTwoAtoms.length; i++) {
                  moleculeTwoAtoms[i].move(-moleculeTwoCOM[0], -moleculeTwoCOM[1], -moleculeTwoCOM[2])
              }
              // Get coordinates of b
              double[] bCoords = b.getXYZ().copy().get()
              // Get rotation matrix to align z-axis with dipole moment (same as above so that the dipoles interact)
              rotation = getRotationBetween(bCoords, zAxis)
              // Get a reference to the positions of all moleculeTwoAtoms
              double[] moleculeTwoAtomPositions = new double[moleculeTwoAtoms.length * 3]
              for(int i = 0; i < moleculeTwoAtoms.length; i++){
                  moleculeTwoAtomPositions[i*3] = moleculeTwoAtoms[i].getX()
                  moleculeTwoAtomPositions[i*3 + 1] = moleculeTwoAtoms[i].getY()
                  moleculeTwoAtomPositions[i*3 + 2] = moleculeTwoAtoms[i].getZ()
              }
              applyRotation(moleculeTwoAtomPositions, rotation)
              // Move atoms into rotated positions
              for(int i = 0; i < moleculeTwoAtoms.length; i++){
                  moleculeTwoAtoms[i].setXYZ(moleculeTwoAtomPositions[i*3],
                          moleculeTwoAtomPositions[i*3 + 1],
                          moleculeTwoAtomPositions[i*3 + 2])
              }
              double[] hBondVector = new double[]{0, 0, a.getZ() - b.getZ() + hBondDist}
              for(int i = 0; i < moleculeTwoAtoms.length; i++){
                  moleculeTwoAtoms[i].move(hBondVector)
              }
              // Minimize the energy of the system subject to a harmonic restraint on the distance between the two atoms
              try {
                  forceFieldEnergy.getCoordinates(x)
                  double e = forceFieldEnergy.energy(x, false)
                  if (e > 100000){
                      throw Exception as Throwable
                  }
                  // Set up restraintBond
                  BondType restraint = new BondType(new int[]{a.getAtomicNumber(), b.getAtomicNumber()},
                          1000.0,
                          hBondDist,
                          BondType.BondFunction.FLAT_BOTTOM_QUARTIC,
                          flatBottomRadius)
                  RestraintBond restraintBond = new RestraintBond(a, b,
                          null,
                          false,
                          0.0, 0.0,
                          null).setBondType(restraint)

                  //Minimize DIRECT and MUTUAL
                  forceFieldEnergy.getPmeNode().setPolarization(Polarization.DIRECT)
                  Minimize minEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
                  minEngine.minimize(1.0, maxIter)
                  forceFieldEnergy.getPmeNode().setPolarization(Polarization.MUTUAL)
                  minEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
                  minEngine.minimize(eps, maxIter)

                  // Delete restraintBond
                  a.getBonds().remove(restraintBond)
                  b.getBonds().remove(restraintBond)
                  a.update()
                  b.update()
                  activeAssembly.bondList.remove(restraintBond)
                  activeAssembly.update()
                  // Store and log the minimized energy
                  forceFieldEnergy.getCoordinates(x)
                  e = forceFieldEnergy.energy(x, false) - monomerEnergy
                  energies.add(e)
                  logger.info(" Energy of trial " + (loopCounter-1) + ": " + e)
                  // Superpose should be run after this to get N^2 rmsds and
                  // then similar structures can be pruned
                  lowestEnergyQueue.add(new StateContainer(new AssemblyState(activeAssembly), e))
              } catch (Exception e) {
                  logger.warning(" Minimization failed. No state will be saved.")
                  //e.printStackTrace()
              }
          }
          zAxis[2] = 1 // Reset z-axis for mol 1 alignment
      }

      logger.info("\n ------------------------- End of Trials -------------------------")
      logger.info(" Lowest energy configuration: " + energies.min())

      // Set up a system XYZ file filter
      File saveLocation = new File(FilenameUtils.removeExtension(filename) + ".arc")
      logger.info(" Logging structures into: " + saveLocation)
      XYZFilter xyzFilter = new XYZFilter(saveLocation,
              activeAssembly,
              activeAssembly.getForceField(),
              activeAssembly.properties)

      // Write to file depending on options
      if(saveAll){
          logger.info(" Saving all " + lowestEnergyQueue.size() + " conformations.")
          int count = 0
          while(!lowestEnergyQueue.empty){ // For loop didn't work?
              StateContainer stateContainer = lowestEnergyQueue.removeFirst()
              AssemblyState state = stateContainer.getState()
              double e = stateContainer.energy
              logger.info(" Writing to file. Configuration #" + (count+1) + " energy: " + e)
              count++
              state.revertState()
              xyzFilter.writeFile(saveLocation, true)
          }
      } else if (saveNumStates != -1 && saveNumStates > 0) {
          int temp = saveNumStates
          saveNumStates = Math.min(saveNumStates, lowestEnergyQueue.size())
          if(temp != saveNumStates) {
              logger.warning(" Input --saveNumStates is greater than number of interactions. Saving all " +
                      "interactions")
          } else {
              logger.info(" Saving out lowest " + saveNumStates + " states.")
          }
          logger.info(" \n")

          for(int i = 0; i < saveNumStates; i++){
              StateContainer stateContainer = lowestEnergyQueue.removeFirst()
              AssemblyState state = stateContainer.getState()
              double e = stateContainer.energy
              logger.info(" Writing to file. Configuration #" + (i+1) + " energy: " + e)
              state.revertState()
              xyzFilter.writeFile(saveLocation, true)
          }
      } else {
          logger.info(" Saving conformations within " + saveEnergyCutoff + " kcal/mol of lowest energy conformation")
          int counter = 1
          double eInit = 10000000
          do {
              StateContainer stateContainer = lowestEnergyQueue.removeFirst()
              AssemblyState state = stateContainer.getState()
              double eTemp = stateContainer.energy
              eInit = counter == 1 ? eTemp : eInit
              if(Math.abs(eTemp - eInit) <= saveEnergyCutoff) {
                  logger.info(" Writing to file. Configuration #" + counter + " energy: " + eTemp)
                  state.revertState()
                  xyzFilter.writeFile(saveLocation, true)
              } else {
                  break
              }
              counter++
          } while(!lowestEnergyQueue.empty)
      }
      return this
  }

    static double[] getCOM(Molecule m){
        // Get center of mass of moleculeOneAtoms
        double[] moleculeOneCOM = new double[3]
        double totalMass = 0.0
        for(Atom s: m.getAtomList()){
            double[] pos = s.getXYZ().get()
            moleculeOneCOM[0] += pos[0] * s.getMass()
            moleculeOneCOM[1] += pos[1] * s.getMass()
            moleculeOneCOM[2] += pos[2] * s.getMass()
            totalMass += s.getMass()
        }
        totalMass = 1 / totalMass
        moleculeOneCOM[0] *= totalMass
        moleculeOneCOM[1] *= totalMass
        moleculeOneCOM[2] *= totalMass

        return moleculeOneCOM
    }

    static double[][] getRotationBetween(double[] v1, double[] v2){
        // Normalize v1 and v2
        double v1Norm = 1/ Math.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2])
        double v2Norm = 1 / Math.sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2])
        for(int i = 0; i < 3; i++){
            v1[i] *= v1Norm
            v2[i] *= v2Norm
        }
        // Cross product between targetVector and z-axis
        double[] crossProduct = new double[3]
        crossProduct[0] = v1[1] * v2[2] - v1[2] * v2[1]
        crossProduct[1] = v1[2] * v2[0] - v1[0] * v2[2]
        crossProduct[2] = v1[0] * v2[1] - v1[1] * v2[0]
        // Normalize cross product
        double crossProductNorm = 1 / Math.sqrt(crossProduct[0] * crossProduct[0] + crossProduct[1] * crossProduct[1] + crossProduct[2] * crossProduct[2])
        for(int i = 0; i < 3; i++){
            crossProduct[i] *= crossProductNorm
        }
        // Dot product between v1 and z-axis
        double dotProduct = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
        // Angle between v1 and z-axis
        double theta = Math.acos(dotProduct)
        double[] u = crossProduct
        double[][] rotation = new double[3][3];
        rotation[0][0] = cos(theta) + (u[0] * u[0]) * (1-cos(theta));
        rotation[0][1] = u[0]*u[1] * (1-cos(theta)) - u[2]*sin(theta);
        rotation[0][2] = u[0]*u[2] * (1-cos(theta)) + u[1]*sin(theta);
        rotation[1][0] = u[1]*u[0] * (1-cos(theta)) + u[2]*sin(theta);
        rotation[1][1] = cos(theta) + (u[1] * u[1]) * (1-cos(theta));
        rotation[1][2] = u[1]*u[2] * (1-cos(theta)) - u[0]*sin(theta);
        rotation[2][0] = u[2]*u[0] * (1-cos(theta)) - u[1]*sin(theta);
        rotation[2][1] = u[2]*u[1] * (1-cos(theta)) + u[0]*sin(theta);
        rotation[2][2] = cos(theta) + (u[2] * u[2]) * (1-cos(theta));
        return rotation
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

