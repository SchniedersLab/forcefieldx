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
import edu.rit.pj.Comm
import ffx.algorithms.AlgorithmListener
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.optimize.Minimize
import ffx.numerics.math.HilbertCurveTransforms
import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.RestraintBond
import ffx.potential.nonbonded.pme.Polarization
import ffx.potential.parameters.BondType
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Molecule
import ffx.potential.bonded.Bond
import ffx.algorithms.optimize.TorsionSearch

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
            description = 'Gradient cutoff for minimization.')
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
            description = 'Initial h-bond distance in angstroms.')
    double hBondDist = 2.0

    /**
     * --flatBottomRadius
     */
    @Option(names = ['--flatBottomRadius', "--fbr"], paramLabel = '0.5',
            description = 'Radius of flat bottom bond restraint potential in angstroms.')
    double flatBottomRadius = 0.5

    /**
     * --saveNumStates
     */
    @Option(names = ['--saveNumStates', "--sns"], paramLabel = '-1',
            description = 'Save this many of the lowest energy states. This is not the default mode.')
    int saveNumStates = -1

    /**
     * --saveEnergyCutoff
     */
    @Option(names = ['--saveEnergyCutoff', "--sec"], paramLabel = '0.0',
            description = 'Cutoff conformations that have an energy (in kcal/mol) less then this cutoff. This is not the default mode.')
    double saveEnergyCutoff = 0.0

    /**
     * --priorTorsionScan
     */
    @Option(names = ['--priorTorsionScan', "--pts"], paramLabel = 'false', defaultValue = 'false',
            description = 'After minimization, statically scan torsions after minimization to find the lowest energy conformation.')
    private boolean priorTorsionScan = false

    /**
     * --intermediateTorsionScan
     */
    @Option(names = ['--intermediateTorsionScan', "--its"], paramLabel = 'false', defaultValue = 'false',
            description = 'During sampling, statically scan torsions after direct minimization to find the lowest energy conformation.')
    private boolean intermediateTorsionScan = false

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

      // Parallelism
      Comm world = Comm.world()

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
          || a.getAtomType().atomicNumber == 16 || a.getAtomType().atomicNumber == 17){ // N,O,F,P,S,Cl
              if(a.getMoleculeNumber() == activeAssembly.getMoleculeNumbers()[0]) {
                  moleculeOneFiltered.add(a)
              } else{
                  moleculeTwoFiltered.add(a)
              }

              // Searching for bonded h's only if we are excluding H's that aren't bonded to electronegative atoms
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
          else if(a.getAtomType().atomicNumber == 1 && !excludeH){ // Add all H's
              if(a.getMoleculeNumber() == activeAssembly.getMoleculeNumbers()[0]) {
                  moleculeOneFiltered.add(a)
              } else{
                  moleculeTwoFiltered.add(a)
              }
          }
      }

      // Make worker assignment adjustments based on moleculeOneFiltered atoms
      if(world.size() > 1){
          // Log total number of atoms
          logger.info("\n Total number of atoms: " + moleculeOneFiltered.size())
          double[] myStates = new double[moleculeOneFiltered.size()]
          for(int i = 0; i < moleculeOneFiltered.size(); i++)
          {
              int assignedWorker = i % world.size()
              if(world.rank() == assignedWorker){
                  myStates[i] = 1 // one-hot encoding assignments of atoms to workers
              }
          }
          // Remove atoms that are not assigned to this worker
          for(int i = moleculeOneFiltered.size() - 1; i >= 0; i--){
              if(myStates[i] == 0){
                  moleculeOneFiltered.remove(i)
              }
          }
          // Log each of the atoms assigned to this worker
          logger.info(" Worker " + world.rank() + " assigned " + moleculeOneFiltered.size() + " atoms.")
          for(Atom a: moleculeOneFiltered){
              logger.info(" " + a.getAtomType().toString())
          }
          // Return if no atoms are assigned to this worker
          if(moleculeOneFiltered.size() == 0){
              logger.info(" Worker " + world.rank() + " has no atoms assigned to it. Exiting.")
              return this
          }
      }

      // Set up forceFieldEnergy
      ForceFieldEnergy forceFieldEnergy = activeAssembly.getPotentialEnergy()
      double[] x = new double[forceFieldEnergy.getNumberOfVariables()]
      forceFieldEnergy.getCoordinates(x)

      // Minimize Monomers, this step is very important since it affects all energies
      double[] minE = minimizeEachSetOfMolecules(activeAssembly, forceFieldEnergy, x, algorithmListener,
              moleculeOneAtoms as ArrayList<Atom>, moleculeTwoAtoms as ArrayList<Atom>, eps, maxIter, priorTorsionScan)
      double monomerEnergy = minE[0]
      double monomerEnergy2 = minE[1]
      double dimerEnergy = minE[2]

      // Log potentials
      logger.info(String.format("\n %-29s%12.7f kcal/mol", "Monomer energy 1:", monomerEnergy))
      monomerEnergy += monomerEnergy2
      logger.info(String.format(" %-29s%12.7f kcal/mol", "Monomer energy 2:", monomerEnergy2))
      logger.info(String.format(" %-29s%12.7f kcal/mol", "Dimer energy:", dimerEnergy))
      logger.info(String.format(" %-29s%12.7f kcal/mol", "Initial Binding Energy:", dimerEnergy - monomerEnergy))

      // Making the MinMax priority queue that will expel the largest entry when it reaches its maximum size
      MinMaxPriorityQueue<StateContainer> lowestEnergyQueue = MinMaxPriorityQueue.
              maximumSize(moleculeOneFiltered.size() * moleculeTwoFiltered.size()).create()
      if(saveNumStates == -1 && saveEnergyCutoff == 0.0) {
          logger.info("\n Saving all conformations.")
      }
      else if (saveNumStates != -1 && saveNumStates > 0 && saveEnergyCutoff == 0.0) {
          lowestEnergyQueue = MinMaxPriorityQueue.maximumSize(saveNumStates).create()
          logger.info("\n Saving the " + saveNumStates + " lowest energy conformations.")
      }
      else {
          saveEnergyCutoff = Math.abs(saveEnergyCutoff)
          logger.info("\n Saving conformations within " + saveEnergyCutoff + " kcal/mol of lowest energy conformation")
      }

      cycleConfigurations(activeAssembly, x, forceFieldEnergy, algorithmListener, molecules,
       moleculeOneFiltered, moleculeOneAtoms, // Parallel workers depend on the ordering
       moleculeTwoFiltered, moleculeTwoAtoms, // of one and two because of loop ordering
       lowestEnergyQueue, monomerEnergy)

      // Set up a system XYZ file filter
      String extension = ".arc"
      if(world.size() > 1){
          extension = "_rank" + world.rank() + extension
      }
      File saveLocation = new File(FilenameUtils.removeExtension(filename) + extension)
      logger.info(" Logging structures into: " + saveLocation)
      XYZFilter xyzFilter = new XYZFilter(saveLocation,
              activeAssembly,
              activeAssembly.getForceField(),
              activeAssembly.properties)

      // Write to file depending on options
      if(saveNumStates == -1 && saveEnergyCutoff == 0.0){
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
      }
      else if (saveNumStates != -1 && saveNumStates > 0 && saveEnergyCutoff == 0.0) {
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
      }
      else {
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

    /**
     * Gets the center of mass of a molecule
     * @param m
     * @return x,y,z coordinates of center of mass
     */
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

    /**
     * Gets the rotation matrix between two vectors
     * @param v1
     * @param v2
     * @return rotation matrix
     */
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
        // Define quaternion from axis-angle
        double[] quaternion = new double[4];
        quaternion[0] = cos(theta/2);
        quaternion[1] = u[0] * sin(theta/2);
        quaternion[2] = u[1] * sin(theta/2);
        quaternion[3] = u[2] * sin(theta/2);
        // Normalize quaternion
        double quaternionNorm = 1 / Math.sqrt(quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1]
                + quaternion[2] * quaternion[2] + quaternion[3] * quaternion[3]);
        for(int i = 0; i < 4; i++){ quaternion[i] *= quaternionNorm; }
        // Useful storage
        double q1q1 = quaternion[1] * quaternion[1];
        double q2q2 = quaternion[2] * quaternion[2];
        double q3q3 = quaternion[3] * quaternion[3];
        double q0q1 = quaternion[0] * quaternion[1];
        double q0q2 = quaternion[0] * quaternion[2];
        double q0q3 = quaternion[0] * quaternion[3];
        double q1q2 = quaternion[1] * quaternion[2];
        double q1q3 = quaternion[1] * quaternion[3];
        double q2q3 = quaternion[2] * quaternion[3];
        // Quaternion rotation matrix
        double[][] rotation = new double[3][3];
        rotation[0][0] = 1 - 2 * (q2q2 + q3q3);
        rotation[0][1] = 2 * (q1q2 - q0q3);
        rotation[0][2] = 2 * (q1q3 + q0q2);
        rotation[1][0] = 2 * (q1q2 + q0q3);
        rotation[1][1] = 1 - 2 * (q1q1 + q3q3);
        rotation[1][2] = 2 * (q2q3 - q0q1);
        rotation[2][0] = 2 * (q1q3 - q0q2);
        rotation[2][1] = 2 * (q2q3 + q0q1);
        rotation[2][2] = 1 - 2 * (q1q1 + q2q2);

        return rotation
    }

    //TODO: Figure out a good way of dealing with the fact that the GPU does not recognize the switch
    // between polarization types in the CPU code does in the next two methods

    /**
     * Does multiple minimization of the active assembly with each molecule using different
     * force field parameters (NONE,DIRECT,MUTUAL)
     * @param activeAssembly
     * @param forceFieldEnergy
     * @param x
     * @param algorithmListener
     * @param moleculeOneAtoms
     * @param moleculeTwoAtoms
     * @param eps
     * @param maxIter
     * @return
     */
    static double[] minimizeEachSetOfMolecules(MolecularAssembly activeAssembly,
                                                        ForceFieldEnergy forceFieldEnergy,
                                                        double[] x,
                                                        AlgorithmListener algorithmListener,
                                                        ArrayList<Atom> moleculeOneAtoms,
                                                        ArrayList<Atom> moleculeTwoAtoms,
                                                        double eps, int maxIter, boolean priorTScan) {
        // Monomer one energy
        for(Atom a: moleculeTwoAtoms){ a.setUse(false) }
        Minimize monomerMinEngine
        logger.info("\n --------- Minimize Monomer 1 --------- ")
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE)
        monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
        monomerMinEngine.minimize(1.0, maxIter).getCoordinates(x)
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.DIRECT)
        monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
        monomerMinEngine.minimize(1.0, maxIter).getCoordinates(x)
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.MUTUAL)
        monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
        monomerMinEngine.minimize(eps, maxIter).getCoordinates(x)
        if(priorTScan) { // Don't want to feel other molecule's effects
            logger.info("\n --------- Monomer 1 Static Torsion Scan --------- ")
            TorsionSearch m1TorsionSearch = new TorsionSearch(activeAssembly, activeAssembly.getMoleculeArray()[0], 16, 1)
            m1TorsionSearch.staticAnalysis(0, 100)
            AssemblyState minState = m1TorsionSearch.getStates().get(0)
            minState.revertState();
        }
        for(Atom a: moleculeTwoAtoms){ a.setUse(true) }

        logger.info("\n --------- Monomer 1 Energy Breakdown --------- ")
        for(Atom a: moleculeTwoAtoms){ a.setUse(false) }
        double monomerEnergy = forceFieldEnergy.energy(x, true)
        for(Atom a: moleculeTwoAtoms){ a.setUse(true) }

        // Monomer two energy
        for(Atom a: moleculeOneAtoms){ a.setUse(false) }
        logger.info("\n --------- Minimize Monomer 2 --------- ")
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE)
        monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
        monomerMinEngine.minimize(1.0, maxIter).getCoordinates(x)
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.DIRECT)
        monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
        monomerMinEngine.minimize(1.0, maxIter).getCoordinates(x)
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.MUTUAL)
        monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
        monomerMinEngine.minimize(eps, maxIter).getCoordinates(x)
        if(priorTScan) {
            logger.info("\n --------- Monomer 2 Static Torsion Scan --------- ")
            TorsionSearch m2TorsionSearch = new TorsionSearch(activeAssembly, activeAssembly.getMoleculeArray()[1], 16, 1)
            m2TorsionSearch.staticAnalysis(0, 100)
            AssemblyState minState = m2TorsionSearch.getStates().get(0)
            minState.revertState();
        }
        for(Atom a: moleculeOneAtoms){ a.setUse(true) }

        logger.info("\n --------- Monomer 2 Energy Breakdown --------- ")
        for(Atom a: moleculeOneAtoms){ a.setUse(false) }
        double monomerEnergy2 = forceFieldEnergy.energy(x, true)
        for(Atom a: moleculeOneAtoms){ a.setUse(true) }

        // Minimize Dimer Structure so energy runs on it
        logger.info("\n --------- Minimize Dimer --------- ")
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE)
        monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
        monomerMinEngine.minimize(1.0, maxIter).getCoordinates(x)
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.DIRECT)
        monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
        monomerMinEngine.minimize(1.0, maxIter).getCoordinates(x)
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.MUTUAL)
        monomerMinEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
        monomerMinEngine.minimize(1.0, maxIter).getCoordinates(x)

        // Get energies logged for init structure
        logger.info("\n --------- Complex Structure Energy Breakdown --------- ")
        double dimerEnergy = forceFieldEnergy.energy(x, true)
        double[] energies = new double[]{monomerEnergy, monomerEnergy2, dimerEnergy}
        return energies
    }


    /**
     * Cycles through all possible configurations of the two molecules and minimizes each one
     * with different force field parameters (NONE,DIRECT,MUTUAL). The lowest energy configuration
     * after minimization is saved to a queue of states (lowest energy queue).
     * @param activeAssembly
     * @param x
     * @param forceFieldEnergy
     * @param algorithmListener
     * @param molecules
     * @param moleculeOneFiltered
     * @param moleculeOneAtoms
     * @param moleculeTwoFiltered
     * @param moleculeTwoAtoms
     * @param lowestEnergyQueue
     * @param monomerEnergy
     */
    private void cycleConfigurations(MolecularAssembly activeAssembly,
                                    double[] x,
                                    ForceFieldEnergy forceFieldEnergy,
                                    AlgorithmListener algorithmListener,
                                    Molecule[] molecules,
                                    ArrayList<Atom> moleculeOneFiltered,
                                    Atom[] moleculeOneAtoms,
                                    ArrayList<Atom> moleculeTwoFiltered,
                                    Atom[] moleculeTwoAtoms,
                                    MinMaxPriorityQueue<StateContainer> lowestEnergyQueue,
                                    double monomerEnergy)
    {
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
                logger.info("\n ----- Trial " + loopCounter + " out of " +
                        (moleculeTwoFiltered.size() * moleculeOneFiltered.size()) + " -----")
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
                    forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE)
                    if(intermediateTorsionScan){ // Molecules feel each other
                        forceFieldEnergy.getCoordinates(x)
                        double energyBefore = forceFieldEnergy.energy(x, false)
                        logger.info(" Energy before Tscan: " + energyBefore)

                        logger.info("\n --------- Monomer 1 Static Torsion Scan --------- ")
                        TorsionSearch m1TorsionSearch = new TorsionSearch(activeAssembly, activeAssembly.getMoleculeArray()[0], 16, 1)
                        m1TorsionSearch.staticAnalysis(0, 100)
                        AssemblyState minState = m1TorsionSearch.getStates().get(0)
                        minState.revertState();

                        logger.info("\n --------- Monomer 2 Static Torsion Scan --------- ")
                        TorsionSearch m2TorsionSearch = new TorsionSearch(activeAssembly, activeAssembly.getMoleculeArray()[1], 16, 1)
                        m2TorsionSearch.staticAnalysis(0, 100)
                        minState = m2TorsionSearch.getStates().get(0)
                        minState.revertState();
                    }
                    forceFieldEnergy.getCoordinates(x)
                    double e = forceFieldEnergy.energy(x, false)
                    if (e > 100000){
                        throw Exception as Throwable
                    }
                    // Set up restraintBond
                    BondType restraint = new BondType(new int[]{a.getAtomicNumber(), b.getAtomicNumber()},
                            1000.0,
                            this.hBondDist,
                            BondType.BondFunction.FLAT_BOTTOM_QUARTIC,
                            this.flatBottomRadius)
                    RestraintBond restraintBond = new RestraintBond(a, b,
                            null,
                            false,
                            0.0, 0.0,
                            null).setBondType(restraint)
                    //Minimize NONE and DIRECT and MUTUAL
                    forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE)
                    Minimize minEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
                    minEngine.minimize(1.0, this.maxIter)
                    forceFieldEnergy.getPmeNode().setPolarization(Polarization.DIRECT)
                    minEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
                    minEngine.minimize(1.0, this.maxIter)
                    forceFieldEnergy.getPmeNode().setPolarization(Polarization.MUTUAL)
                    minEngine = new Minimize(activeAssembly, forceFieldEnergy, algorithmListener)
                    minEngine.minimize(this.eps, this.maxIter)
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
                    logger.info(" Binding energy of trial " + (loopCounter-1) + ": " + e)
                    // Superpose should be run after this to get N^2 rmsds and
                    // then similar structures can be pruned
                    lowestEnergyQueue.add(new StateContainer(new AssemblyState(activeAssembly), e))
                } catch (Exception ignored) {
                    logger.warning(" Minimization failed. No state will be saved.")
                    //e.printStackTrace()
                }
            }
            zAxis[2] = 1 // Reset z-axis for mol 1 alignment
        }
        logger.info("\n ------------------------- End of Trials -------------------------")
    }

    /**
     * Implements StateContainer to store the coordinates of a state and its energy
     */
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

