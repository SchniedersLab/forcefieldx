//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.potential.groovy.test

import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.MSNode
import ffx.potential.cli.PotentialScript
import ffx.potential.utils.PACCOMFunctions
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.Superpose
import ffx.utilities.DoubleIndexPair
import org.apache.commons.io.FilenameUtils
import org.apache.commons.lang.ArrayUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.logging.Level
import java.util.stream.Collectors

import static java.lang.String.format

/**
 * The PACCOM script calculates intracrystollographic distances.
 *
 * @author Okimasa OKADA
 * created by Okimasa OKADA 2017/3/31
 * revised by Okimasa OKADA 2019/2/25
 * ported to FFX by Aaron Nessler, Kaleb Bierstedt, and Micheal Schnieders 2020
 * <br>
 * Usage:
 * <br>
 * ffxc test.PACCOM &lt;filename&gt &lt;filename&gt;
 */
@Command(description = " Compare crystal packings based on intermolecular distances.", name = "ffxc test.PACCOM")
class PACCOM extends PotentialScript {

  /**
   * --nm or --numberMolecules Number of molecules to include from each crystal in RMSD comparison.
   */
  @Option(names = ['--nm', '--numberMolecules'], paramLabel = '50',
      description = 'Determines crystal sphere size for comparison.')
  int nMolecules = 50

  /**
   * --ls or --largeSphere Number of molecules in the large sphere to compare with first crystal.
   */
  @Option(names = ['--ls', '--largeSphere'], paramLabel = '500',
      description = 'Determines number of molecules in initial sphere.')
  int largeSphere = 500

  /**
   * --all or --allVsAll Compute RMSD of between each file.
   */
  @Option(names = ['--all', '--allVsAll'], paramLabel = "false", defaultValue = "false",
      description = 'Calculate the RMSD between each sphere.')
  private boolean allVsAll = false

  /**
   * --op or --originalPACCOM Use Algorithm supplied by original PACCOM.
   */
  @Option(names = ['--op', '--originalPACCOM'], paramLabel = "false", defaultValue = "false",
      description = 'Utilize a version of PACCOM that mirrors original.')
  private static boolean original = false

  /**
   * -s or --save Save out individual XYZ/PDB files for each created sphere.
   */
  @Option(names = ['-s', '--save'], paramLabel = "false", defaultValue = "false",
      description = 'Save each sphere as an aperiodic system to a PDB file.')
  private boolean saveFiles = false

  // TODO: implement closestDistance
  /**
   * --cd or --closestDistance Neighbor molecules calculated via closest atom (not currently implemented).
   */
//    @Option(names = ['--cd', '--closestDistance'], paramLabel = "false", defaultValue = "false",
//            description = 'Neighbors calculated via closest atom.')
//    private static boolean closestDistance = false

  /**
   * --nh or --noHydrogens Perform comparison without hydrogens.
   */
  @Option(names = ['--nh', '--noHydrogens'], paramLabel = "false", defaultValue = "false",
      description = 'Crystal RMSD calculated without hydrogen atoms.')
  private static boolean noHydrogens = false

  /**
   * Select individual atoms to match and calculate RMSD rather than using all atoms.
   */
  @Option(names = ['--ca', '--centerAtoms'], arity = "1..*", paramLabel = "atom integers",
      description = 'Specify atoms to match and calculate RMSD. Otherwise, all atom.')
  private int[] centerAtoms = null

  /**
   * The final argument(s) should be two or more filenames.
   */
  @Parameters(arity = "2..*", paramLabel = "files",
      description = 'Atomic coordinate files to compare in XYZ format.')
  List<String> filenames = null

  private MolecularAssembly[] assemblies

  public double[][] coordinates = null

  private File baseDir = null

  /**
   * JUnit Testing Variables
   */
  public double[][] rmsdValues = null

  /**
   * PACCOM Constructor.
   */
  PACCOM() {
    this(new Binding())
  }

  /**
   * PACCOM Constructor.
   * @param binding Groovy Binding to use.
   */
  PACCOM(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  PACCOM run() {
    System.setProperty("vdwterm", "false")
    // Read in structures
    if (!init()) {
      return
    }
    // Ensure file exists/create active molecular assembly
    assemblies = new MolecularAssembly[filenames.size()]
    if (filenames != null && filenames.size() > 0) {
      for (int i = 0; i < filenames.size(); i++) {
        assemblies[i] = potentialFunctions.open(filenames.get(i))
      }
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return
    } else {
      assemblies = [activeAssembly]
    }

    Binding binding = getBinding()
    baseDir = assemblies[0].getFile()
    potentialFunctions = (PotentialsFunctions) binding.getVariable("functions")

    int numAssemblies = assemblies.size()
    // TODO: Make more robust? This assumes atom orders are the same in all molecules... (noHydrogens)
    List<Atom> exampleAtoms = assemblies[0].getMolecules().get(0).getAtomList()
    int nAtoms = exampleAtoms.size()
    int numHydrogens = 0
    MolecularAssembly[] expandedAssemblies = new MolecularAssembly[numAssemblies]
    MolecularAssembly[] inflatedAssemblies = new MolecularAssembly[numAssemblies]
    // If comparison atoms are not specified use all atoms from the molecule.
    int[] comparisonAtoms
    if (centerAtoms != null) {
      comparisonAtoms = centerAtoms.clone()
      for (int i = 0; i < comparisonAtoms.size(); i++) {
        comparisonAtoms[i] = comparisonAtoms[i] - 1
        if (comparisonAtoms[i] < 0 ||
            comparisonAtoms[i] > assemblies[0].getMolecules().get(0).getAtomList(false).size()) {
          logger.severe(" Selected atoms are outside of molecular size.")
        }
      }
    } else {
      if (noHydrogens) {
        int[] hydrogenIndices = 0
        numHydrogens = 0
        for (int i = 0; i < exampleAtoms.size(); i++) {
          if (exampleAtoms.get(i).isHydrogen()) {
            hydrogenIndices[numHydrogens++] = i
          }
        }
        comparisonAtoms = new int[exampleAtoms.size() - numHydrogens]
        for (int i = 0; i < comparisonAtoms.size(); i++) {
            if (!ArrayUtils.contains(hydrogenIndices, i)) {
                comparisonAtoms[i] = i
            }
        }

      } else {
        comparisonAtoms = new int[exampleAtoms.size()]
        for (int i = 0; i < comparisonAtoms.size(); i++) {
          comparisonAtoms[i] = i
        }
      }
    }

    int compareAtomsSize = comparisonAtoms.size()
    logger.info(format(" Number of atoms to compare: %d", compareAtomsSize))
    for (int integerVal : comparisonAtoms) {
      logger.fine(format(" Comparing atom #: %d", integerVal + 1))
    }

    // Values for Printing RMSD at the end.
    rmsdValues = new double[numAssemblies][numAssemblies]
    String row = new String()

    int massIndex
    int coordIndex
//        int mnum = 1 // TODO: determine better value than nMolecules
    int numPDBs = 0

    double[] targetCoords
    double[][] assemblyCoords
    //for (int m = 0; m < mnum; m++) {
    // Loop over each crystal, move to origin, orient to x and x-y
//        if (m == 0) {

    DoubleIndexPair[][] expandedPairs = new DoubleIndexPair[numAssemblies][largeSphere]

    for (int i = 0; i < numAssemblies; i++) {
      String fileName = FilenameUtils.getBaseName(assemblies[i].getFile().getName())
      inflatedAssemblies[i] =
          PACCOMFunctions.generateBaseSphere(assemblies[i], largeSphere, expandedPairs[i], true)
      // TODO replace boolean with original
      MolecularAssembly currentAssembly = inflatedAssemblies[i]
      logger.info(format(" Number of molecules in currentAssembly: %d",
          currentAssembly.getMolecules().size()))
      if (saveFiles && logger.isLoggable(Level.FINE)) {
        PACCOMFunctions.saveAssemblyPDB(currentAssembly, fileName + "_inflated", numPDBs++)
      }

    }
//        } else {
//            for (int i = 1; i < numAssemblies; i++) {
//                String fileName = FilenameUtils.getBaseName(assemblies[i].getFile().getName())
//
//                inflatedAssemblies[i] = PACCOMFunctions.generateBaseSphere(assemblies[i], largeSphere, true)
//                // TODO replace boolean with original
//                MolecularAssembly currentAssembly = inflatedAssemblies[i]
//                for (int l = 0; l < currentAssembly.getMolecules().size(); l++) {
//                    //logger.info(format("Molecule %d" + currentAssembly.getMolecules().get(l).getAtomList().get(comparisonAtoms[0]), l))
//                }
//                logger.info(format(" Number of molecules in currentAssembly: %d", currentAssembly.getMolecules().size()))
//                if (saveFiles && logger.isLoggable(Level.FINE)) {
//                    PACCOMFunctions.saveAssemblyPDB(currentAssembly, fileName +format("_%d_%d", m, i) + "inflated", numPDBs++)
//                }
//                minIndices[i] = PACCOMFunctions.moleculeIndicesNear(currentAssembly, m, i, minIndices, nMolecules, original)
//            }
//        }

    //**            Move to origin         **//
    ArrayList<Atom> list1Center
    Atom atomN11
    Atom atomN21
    Atom atomN31
    assemblyCoords = new double[numAssemblies][inflatedAssemblies[0].getAtomList().size() * 3]
    for (int i = 0; i < numAssemblies; i++) {
      String fileName = FilenameUtils.getBaseName(assemblies[i].getFile().getName())
      MolecularAssembly currentAssembly = inflatedAssemblies[i]
      int indexOfClosestMolecule = expandedPairs[i][0].getIndex()
      MSNode centerMolecule = currentAssembly.getMolecules().get(indexOfClosestMolecule)
      ArrayList<Atom> centerMolAtomList = centerMolecule.getAtomList()
      ArrayList<Atom> assemblyAtomList = currentAssembly.getAtomList()
      int assemblyAtomsSize = assemblyAtomList.size()
      double[] assemblyMasses = new double[assemblyAtomsSize]
      // Record coordinates for all atoms within the assembly in a manner the Superpose class will handle.
      massIndex = 0
      coordIndex = 0
      for (Atom atom : assemblyAtomList) {
        assemblyMasses[massIndex++] = atom.getMass()
        assemblyCoords[i][coordIndex++] = atom.getX()
        assemblyCoords[i][coordIndex++] = atom.getY()
        assemblyCoords[i][coordIndex++] = atom.getZ()
      }

      // Begin translation of center molecule to origin.

      double[] coordsToMove = new double[3]
      double[] centerMasses = new double[1]
      centerMasses[0] = 1.0

      Atom atomN1 = centerMolAtomList.get(comparisonAtoms[0])
      coordsToMove[0] = atomN1.getX()
      coordsToMove[1] = atomN1.getY()
      coordsToMove[2] = atomN1.getZ()

      logger.fine(format(" Atom coordinates pre translation:"))
      logger.fine(format(" N1: " + atomN1.toString()))
      logger.fine(format(" N2: " + centerMolAtomList.get(comparisonAtoms[1]).toString()))
      logger.fine(format(" N3: " + centerMolAtomList.get(comparisonAtoms[2]).toString()))

      // Determine/apply translation necessary to get CoM for middle molecule to origin
      double[] translation = Superpose.calculateTranslation(coordsToMove, centerMasses)
      String translationString = new String()
      for (double valueT : translation) {
        translationString += format("%16.8f", valueT) + "\t"
      }
      logger.info(format(" Translation Vector: (Move to Origin)\n %s", translationString))
      Superpose.applyTranslation(assemblyCoords[i], translation)

      // Translation to origin finished. Begin first rotation:

      atomN1 = centerMolAtomList.get(comparisonAtoms[0])
      Atom atomN2 = centerMolAtomList.get(comparisonAtoms[1])
      Atom atomN3 = centerMolAtomList.get(comparisonAtoms[2])

      coordsToMove = new double[nMolecules * 3]
      coordsToMove[0] = atomN1.getX()
      coordsToMove[1] = atomN1.getY()
      coordsToMove[2] = atomN1.getZ()
      coordsToMove[3] = atomN2.getX()
      coordsToMove[4] = atomN2.getY()
      coordsToMove[5] = atomN2.getZ()
      coordsToMove[6] = atomN3.getX()
      coordsToMove[7] = atomN3.getY()
      coordsToMove[8] = atomN3.getZ()

      targetCoords = PACCOMFunctions.standardOrientation(atomN1, atomN2, atomN3)

      centerMasses = new double[nMolecules]
      Arrays.fill(centerMasses, 1.0)

      double[][] rotation = Superpose.calculateRotation(targetCoords, coordsToMove, centerMasses)
      // Print out first rotation matrix
      String rotString = new String()
      for (double[] rotRow : rotation) {
        rotString += "\n  "
        for (double rotCol : rotRow) {
          rotString += format("%16.8f", rotCol) + "\t"
        }
      }
      logger.info(" First Rotation Matrix:" + rotString)

      Superpose.applyRotation(assemblyCoords[i], rotation)

      // Update expanded assemblies to include rotation.
      coordIndex = 0
      for (Atom atom : assemblyAtomList) {
        double[] newXYZ = new double[] {assemblyCoords[i][coordIndex++],
            assemblyCoords[i][coordIndex++],
            assemblyCoords[i][coordIndex++]}
        atom.setXYZ(newXYZ)
      }

      if (saveFiles && logger.isLoggable(Level.FINE)) {
        PACCOMFunctions.saveAssemblyPDB(currentAssembly, fileName + "_1stRot", numPDBs++)
      }
      // Check RMSD of center molecules in shell compared to first
      if (i == 0) {
        int index = expandedPairs[i][0].getIndex()
        list1Center = inflatedAssemblies[0].getMolecules().get(index).getAtomList()
        atomN11 = list1Center.get(comparisonAtoms[0])
        atomN21 = list1Center.get(comparisonAtoms[1])
        atomN31 = list1Center.get(comparisonAtoms[2])
      }

      int index = expandedPairs[i][0].getIndex()
      ArrayList<Atom> list2Center = inflatedAssemblies[i].getMolecules().get(index).getAtomList()

      Atom atomN12 = list2Center.get(comparisonAtoms[0])
      Atom atomN22 = list2Center.get(comparisonAtoms[1])
      Atom atomN32 = list2Center.get(comparisonAtoms[2])

      double[][] checkCoords = new double[2][compareAtomsSize * 3]
      checkCoords[0][0] = atomN11.getX()
      checkCoords[0][1] = atomN11.getY()
      checkCoords[0][2] = atomN11.getZ()
      checkCoords[0][3] = atomN21.getX()
      checkCoords[0][4] = atomN21.getY()
      checkCoords[0][5] = atomN21.getZ()
      checkCoords[0][6] = atomN31.getX()
      checkCoords[0][7] = atomN31.getY()
      checkCoords[0][8] = atomN31.getZ()

      checkCoords[1][0] = atomN12.getX()
      checkCoords[1][1] = atomN12.getY()
      checkCoords[1][2] = atomN12.getZ()
      checkCoords[1][3] = atomN22.getX()
      checkCoords[1][4] = atomN22.getY()
      checkCoords[1][5] = atomN22.getZ()
      checkCoords[1][6] = atomN32.getX()
      checkCoords[1][7] = atomN32.getY()
      checkCoords[1][8] = atomN32.getZ()

      double[] checkMass = new double[] {1.0, 1.0, 1.0}

      double tempRMSD = Superpose.rmsd(checkCoords[0], checkCoords[1], checkMass)
      logger.info(format("Post-rotation RMSD between center molecules: %16.8f", tempRMSD))

//            logger.fine(format(" Center N1 atoms for %d shell.", i))
//            for (int l : minIndices[i]) {
//                Atom atom = inflatedAssemblies[i].getMolecules().get(l).getAtomList().get(comparisonAtoms[0])
//                logger.fine(format("HETATM %4d  C   100 A %3d   %8.3f%8.3f %8.3f", l + 1, l + 1, atom.getX(), atom.getY(), atom.getZ()))
//                PACCOMFunctions.saveAssemblyPDB(inflatedAssemblies[0], "list1", numPDBs++)
//            }
    }

    /**     PERFORM nMolecules MOLECULE SUPERPOSITION     **/
    // At this point, all expanded assemblies should have their center molecules aligned at the origin.
    logger.info("\n\n Begin Second Rotation\n\n")

    for (int i = 0; i < numAssemblies; i++) {
      int[][] similarMolecules = new int[numAssemblies][nMolecules]
      // Need to have array initialized with impossible values
      String fileName = FilenameUtils.getBaseName(assemblies[i].getFile().getName())
      MolecularAssembly currentAssembly = inflatedAssemblies[i]
      // Find molecules closest to target assembly in current assembly
//                Arrays.fill(similarMolecules[ii], -1)
//                int targetIndex = 0
//                for (int j = 0; j < nMolecules; j++) {
//                    double[] targetCenter = new double[3]
//                    targetAssembly.getMolecules().get(minIndices[i][j]).getAtomList().get(comparisonAtoms[0]).getXYZ(targetCenter)
//                    int index = -1
//                    double minDist = Double.MAX_VALUE
//                    double value
//                    for (int k = 0; k < currentAssembly.getMolecules().size(); k++) {
//                        boolean duplicate = false
//                        double[] nextCenter = new double[3]
//                        currentAssembly.getMolecules().get(k).getAtomList().get(comparisonAtoms[0]).getXYZ(nextCenter)
//                        value = distance.compute(targetCenter, nextCenter)
//                        if (value < minDist) {
//                            for (int prevKey : similarMolecules[ii]) {
//                                if (prevKey == k) {
//                                    duplicate = true
//                                    break
//                                }
//                            }
//                            logger.fine(format(" Key: " + k + " Value: " + value + " is duplicate? " + duplicate))
//                            if (!duplicate) {
//                                minDist = value
//                                index = k
//                            }
//                        }
//                    }
//                    if (minDist > 1.0) {
//                        logger.warning(format(" Matched molecules are greater than 1.0 Å (%16.8f)", minDist))
//                    }
//                    similarMolecules[ii][targetIndex++] = index
//                }
      // Prepare coordinates for Superpose algorithm
//                Atom atomN1 = assemblyAtomList.get(nAtoms * similarMolecules[ii][0] + comparisonAtoms[0])
//                Atom atomN2 = assemblyAtomList.get(nAtoms * similarMolecules[ii][1] + comparisonAtoms[0])
//                Atom atomN3 = assemblyAtomList.get(nAtoms * similarMolecules[ii][2] + comparisonAtoms[0])

      Atom atomN1 = currentAssembly.getAtomList().
          get(nAtoms * expandedPairs[i][0].getIndex() + comparisonAtoms[0])
      Atom atomN2 = currentAssembly.getAtomList().
          get(nAtoms * expandedPairs[i][1].getIndex() + comparisonAtoms[0])
      Atom atomN3 = currentAssembly.getAtomList().
          get(nAtoms * expandedPairs[i][2].getIndex() + comparisonAtoms[0])

      targetCoords = PACCOMFunctions.standardOrientation(atomN1, atomN2, atomN3)

      double[] coordsToMove = new double[compareAtomsSize * 3]
      coordsToMove[0] = atomN1.getX()
      coordsToMove[1] = atomN1.getY()
      coordsToMove[2] = atomN1.getZ()
      coordsToMove[3] = atomN2.getX()
      coordsToMove[4] = atomN2.getY()
      coordsToMove[5] = atomN2.getZ()
      coordsToMove[6] = atomN3.getX()
      coordsToMove[7] = atomN3.getY()
      coordsToMove[8] = atomN3.getZ()

      double[] masses = new double[compareAtomsSize]
      // Oringial PACCOM does not use masses
      Arrays.fill(masses, 1.0)

      logger.finest(
          format(" target: %d coords: %d masses: %d", targetCoords.length, coordsToMove.length,
              masses.length))
      double[][] rotation = Superpose.calculateRotation(targetCoords, coordsToMove, masses)
      String rotString = new String()
      for (double[] rotRow : rotation) {
        rotString += "\n  "
        for (double rotCol : rotRow) {
          rotString += format("%16.8f", rotCol) + "\t"
        }
      }
      logger.info(" Second Rotation Matrix:" + rotString)
      Superpose.applyRotation(assemblyCoords[i], rotation)

      // Update expanded assemblies to include second rotation.
      coordIndex = 0
      for (Atom atom : currentAssembly.getAtomList()) {
        double[] newXYZ = new double[] {assemblyCoords[i][coordIndex++],
            assemblyCoords[i][coordIndex++],
            assemblyCoords[i][coordIndex++]}
        atom.setXYZ(newXYZ)
      }

      if (saveFiles && logger.isLoggable(Level.FINE)) {
        PACCOMFunctions.saveAssemblyPDB(currentAssembly, fileName + "_2ndRot", numPDBs++)
      }
      // Determine molecules closest to shell one for
//                targetIndex = 0
//                for (int j = 0; j < minIndices[i].length; j++) {
//                    double[] targetCenter = new double[3]
//                    targetAssembly.getMolecules().get(minIndices[i][j]).getAtomList().get(comparisonAtoms[0]).getXYZ(targetCenter)
//                    int index = -1
//                    double minDist = Double.MAX_VALUE
//                    double value
//                    for (int k = 0; k < currentAssembly.getMolecules().size(); k++) {
//                        boolean duplicate = false
//                        double[] nextCenter = new double[3]
//                        currentAssembly.getMolecules().get(k).getAtomList().get(comparisonAtoms[0]).getXYZ(nextCenter)
//                        value = distance.compute(targetCenter, nextCenter)
//                        if (value < minDist) {
//                            for (int prevKey : similarMolecules[ii]) {
//                                if (prevKey == k) {
//                                    duplicate = true
//                                    break
//                                }
//                            }
//                            //logger.fine(format(" Key: " + k + " Value: " + value + " is duplicate? " + duplicate))
//                            if (!duplicate) {
//                                minDist = value
//                                index = k
//                            }
//                        }
//                    }
//                    if (minDist > 1.0) {
//                        //logger.warning(format(" Matched molecules are greater than 1.0 Å (%16.8f)", minDist))
//                    }
//                    //logger.fine(format(" nMol Index: %d \tDistance to Matched: %16.8f", index, minDist))
//                    similarMolecules[ii][targetIndex++] = index
//                }

//                expandedAssemblies[ii] = PACCOMFunctions.cutSubSphere(assemblies[ii], currentAssembly, similarMolecules[ii])
      expandedAssemblies[i] =
          PACCOMFunctions.cutSubSphere(assemblies[i], currentAssembly, expandedPairs[i], nMolecules)
      if (saveFiles) {
        PACCOMFunctions.saveAssemblyPDB(expandedAssemblies[i], fileName + "_final", numPDBs++)
      }
    }
    for (int i = 0; i < numAssemblies; i++) {
      String tempTarget = " "
      row += format(" %s\t", filenames.get(i))
      MolecularAssembly targetAssembly = expandedAssemblies[i]
      List<Atom> targetAtomList = targetAssembly.getAtomList()
      for (Atom atom : targetAtomList) {
        if (noHydrogens) {
          if (!atom.isHydrogen()) {
            tempTarget += atom.getAtomType().name
          }
        } else {
          tempTarget += atom.getAtomType().name
        }
      }
      for (int ii = 0; ii < numAssemblies; ii++) {
        MolecularAssembly currentAssembly = expandedAssemblies[ii]
        // Closest nMolecules have been determined for each sphere
        String tempCurrent = " "
        ArrayList<Atom> assemblyAtomList = currentAssembly.getAtomList()
        for (Atom atom : assemblyAtomList) {
          if (noHydrogens) {
            if (!atom.isHydrogen()) {
              tempCurrent += atom.getAtomType().name
            }
          } else {
            tempCurrent += atom.getAtomType().name
          }
        }
        // String comparison to ensure same atom types are being matched.
        if (tempTarget != tempCurrent) {
          for (int j = 0; j < tempTarget.length(); j++) {
            if (tempTarget.charAt(j) != tempCurrent.charAt(j)) {
              logger.warning(" Atom name (" + tempTarget.charAt(j) +
                  ") from target does not match comparison atom (" + tempCurrent.charAt(j) + ").")
            }
          }
        }
        // Perform RMSD calculation
        double[] targetFinalCoords = new double[targetAssembly.getMolecules().size() * nAtoms * 3]
        coordIndex = 0
        for (MSNode molecule : targetAssembly.getMolecules()) {
          for (Atom atom : molecule.getAtomList()) {
            if (noHydrogens) {
              if (!atom.isHydrogen()) {
                targetFinalCoords[coordIndex++] = atom.getX()
                targetFinalCoords[coordIndex++] = atom.getY()
                targetFinalCoords[coordIndex++] = atom.getZ()
              }
            } else {
              targetFinalCoords[coordIndex++] = atom.getX()
              targetFinalCoords[coordIndex++] = atom.getY()
              targetFinalCoords[coordIndex++] = atom.getZ()
            }
          }
        }

        // Record coordinates for all atoms within the assembly in a manner the Superpose class will handle.
        int numMolRMSD = currentAssembly.getMolecules().size()
        double[] assemblyMasses
        if (noHydrogens) {
          assemblyMasses = new double[numMolRMSD * (nAtoms - numHydrogens)]
        } else {
          assemblyMasses = new double[numMolRMSD * nAtoms]
        }
        Arrays.fill(assemblyMasses, 1.0)
        double[] currentAssemblyCoords = new double[numMolRMSD * nAtoms * 3]
        coordIndex = 0
        massIndex = 0
        for (MSNode molecule : currentAssembly.getMolecules()) {
          for (Atom atom : molecule.getAtomList()) {
            if (noHydrogens) {
              if (!atom.isHydrogen()) {
                assemblyMasses[massIndex++] = atom.getMass()
                currentAssemblyCoords[coordIndex++] = atom.getX()
                currentAssemblyCoords[coordIndex++] = atom.getY()
                currentAssemblyCoords[coordIndex++] = atom.getZ()
              }
            } else {
              assemblyMasses[massIndex++] = atom.getMass()
              currentAssemblyCoords[coordIndex++] = atom.getX()
              currentAssemblyCoords[coordIndex++] = atom.getY()
              currentAssemblyCoords[coordIndex++] = atom.getZ()
            }
          }
        }
        if (original) {
          Arrays.fill(assemblyMasses, 1.0)
        }
        rmsdValues[i][ii] = Superpose.rmsd(targetFinalCoords, currentAssemblyCoords, assemblyMasses)
        row += format("%16.8f\t", rmsdValues[i][ii])
      }
      row += format("\n")
    }
    logger.info(row)
    return this
  }

  @Override
  List<Potential> getPotentials() {
    if (assemblies == null) {
      return new ArrayList<Potential>()
    } else {
      return Arrays.stream(assemblies).
          filter {a -> a != null
          }.
          map {a -> a.getPotentialEnergy()
          }.
          filter {e -> e != null
          }.
          collect(Collectors.toList())
    }
  }
}