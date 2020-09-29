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

import ffx.numerics.math.Double3
import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.stream.IntStream

import static ffx.potential.utils.Superpose.*
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.sqrt

/**
 * The Superpose script superposes molecules in an arc/multiple model pdb file (all versus all or one versus all) or in two pdb/xyz files.
 * <br>
 * Usage:
 * <br>
 * ffxc Superpose [options] &lt;filename&gt;
 */
@Command(description = " Superpose frames of a trajectory file and calculate RMSD.", name = "ffxc Superpose")
class Superpose extends PotentialScript {

  /**
   * --aS or --atomSelection The atom selection [HEAVY (0) / ALL (1) / CALPHA (2) / BACKBONE (3)] for the RMSD calculation (CALPHA chooses N1 or N9 for nucleic acids).
   */
  @Option(names = ['--aS', '--atomSelection'], paramLabel = "0", defaultValue = "0",
      description = 'The atom selection [HEAVY (0) / ALL (1) / CALPHA (2) / BACKBONE (3) ] for the RMSD calculation (CALPHA chooses N1 or N9 for nucleic acids).')
  private String atomSelection = "0"

  /**
   * -A or --allvsAll Frames to be compared within the arc file. Select [true] for all versus all comparison; select [false] for one versus all comparison.
   */
  @Option(names = ['-A', '--allvsAll'], paramLabel = "false", defaultValue = "false",
      description = 'Compare all snapshots versus all others, instead of the first snapshot versus all others.')
  private boolean frameComparison = false

  /**
   * --store or --storeMatrix Store the distance matrix from all versus all RMSD calculation on multiple models.
   */
  @Option(names = ['--store', '--storeMatrix'], paramLabel = "false", defaultValue = "false",
      description = 'Store the distance matrix of all versus all RMSD calculation.')
  private boolean storeMatrix = false

  /**
   * -s or --start Atom number where RMSD calculation of structure will begin.
   */
  @Option(names = ['-s', '--start'], paramLabel = "1", defaultValue = "1",
      description = 'Starting atom to include in the RMSD calculation.')
  private int start = 1

  /**
   * -f or --final Atom number where RMSD calculation of structure will end.
   */
  @Option(names = ['-f', '--final'], paramLabel = "nAtoms",
      description = 'Final atom to include in the RMSD calculation.')
  private int finish = Integer.MAX_VALUE

  /**
   * --dRMSD Calculate the dRMSD in addition to RMSD.
   */
  @Option(names = ['--dRMSD'], paramLabel = "nAtoms",
      description = 'Calculate the dRMSD in addtion to RMSD.')
  private boolean dRMSD = false

  /**
   * -w or --write Write out superposed snapshots.
   */
  @Option(names = ['-w', '--write'], paramLabel = "false", defaultValue = "false",
      description = 'Write out superposed snapshots.')
  private boolean writeSnapshots = false

  /**
   * -v or --verbose Print out RMSD information.
   */
  @Option(names = ['-v', '--verbose'], paramLabel = "true", defaultValue = "true",
      description = 'Write out RMSD information.')
  private boolean verbose = true

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1", paramLabel = "files",
      description = 'The atomic coordinate file in PDB or XYZ format.')
  List<String> filenames = null

  public ForceFieldEnergy forceFieldEnergy = null
  private File outFile
  private XYZFilter outputFilter

  double[][] distMatrix

  /**
   * Superpose Constructor.
   */
  Superpose() {
    this(new Binding())
  }

  /**
   * Superpose Constructor.
   * @param binding Groovy Binding to use.
   */
  Superpose(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  Superpose run() {
    if (!init()) {
      return this
    }

    // Turn off non-bonded terms for efficiency.
    System.setProperty("vdwterm", "false")

    if (filenames != null && filenames.size() > 0) {
      MolecularAssembly[] assemblies = [potentialFunctions.open(filenames.get(0))]
      activeAssembly = assemblies[0]
      if (filenames.size() > 1) {
        MolecularAssembly[] assemblies2 = [potentialFunctions.open(filenames.get(1))]
      }
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // Either the SystemFilter for the first parsed file, or the 2nd parsed file.
    SystemFilter systemFilter = potentialFunctions.getFilter()

    forceFieldEnergy = activeAssembly.getPotentialEnergy()
    Atom[] atoms = activeAssembly.getAtomArray()

    int nVars = forceFieldEnergy.getNumberOfVariables()
    double[] x = new double[nVars]
    forceFieldEnergy.getCoordinates(x)

    if (writeSnapshots) {
      String outFileName =
          activeAssembly.getFile().toString().replaceFirst(~/\.(?:xyz|pdb|arc).*$/, "")
      outFileName = outFileName + "_superposed.arc"
      outFileName = potentialFunctions.versionFile(outFileName)

      outFile = new File(outFileName)
      outputFilter = new XYZFilter(outFile, activeAssembly, activeAssembly.getForceField(),
          activeAssembly.getProperties())
    }

    int distMatrixSize = systemFilter.countNumModels()
    distMatrix = new double[distMatrixSize][distMatrixSize]

    if (systemFilter instanceof PDBFilter || systemFilter instanceof XYZFilter) {
      double[] x2 = new double[nVars]
      double[] mass = new double[nVars / 3]

      int nAtoms = atoms.length
      for (int i = 0; i < nAtoms; i++) {
        mass[i] = atoms[i].getMass()
      }

      if (finish > nAtoms - 1) {
        finish = nAtoms - 1
      }
      if (start < 0 || start > finish) {
        start = 0
      }

      // Note that atoms are indexed from 0 to nAtoms - 1.
      if (verbose) {
        logger.info(format("\n Atoms from %d to %d will be considered.", start, finish))
      }

      // Begin streaming the possible atom indices, filtering out inactive atoms.
      IntStream atomIndexStream =
          IntStream.range(start, finish + 1).filter({int i -> return atoms[i].isActive()
          })

      // String describing the selection type.
      String selectionType = "All Atoms"

      // Switch on what type of atoms to select, filtering as appropriate. Support the old integer indices.
      switch (atomSelection.toUpperCase()) {
        case "HEAVY":
        case "0":
          // Filter only for heavy (non-hydrogen) atoms.
          atomIndexStream = atomIndexStream.filter({int i -> atoms[i].isHeavy()
          })
          selectionType = "Heavy Atoms"
          break
        case "ALL":
        case "1":
          // Unmodified stream; we have just checked for active atoms.
          selectionType = "All Atoms"
          break
        case "ALPHA":
        case "2":
          // Filter only for reference atoms: carbons named CA (protein) or nitrogens named N1 or N9 (nucleic acids).
          atomIndexStream = atomIndexStream.filter({int i ->
            Atom ati = atoms[i]
            String atName = ati.getName().toUpperCase()
            boolean proteinReference = atName == "CA" && ati.getAtomType().atomicNumber == 6
            boolean naReference = (atName == "N1" || atName == "N9") &&
                ati.getAtomType().atomicNumber == 7
            return proteinReference || naReference
          })
          selectionType = "C-Alpha Atoms (or N1/N9 for nucleic acids)"
          break
        case "BACKBONE":
        case "3":
          // Filter for only backbone atoms.
          atomIndexStream = atomIndexStream.filter({int i ->
            Atom ati = atoms[i]
            String atName = ati.getName().toUpperCase()
            boolean caReference = atName.equals("CA") && ati.getAtomType().atomicNumber == 6
            boolean cReference = atName.equals("C") && ati.getAtomType().atomicNumber == 6
            boolean nReference = atName.equals("N") && ati.getAtomType().atomicNumber == 7
            return caReference || cReference || nReference
          })
          selectionType = "C, C-Alpha, and N backbone atoms."
          break
        default:
          logger.severe(format(
              " Could not parse %s as an atom selection! Must be ALL, HEAVY, ALPHA or BACKBONE.",
              atomSelection))
          break
      }

      if (verbose) {
        logger.info(" Superpose selection criteria: " + selectionType + "\n")
      }

      // Indices of atoms used in alignment and RMSD calculations.
      int[] usedIndices = atomIndexStream.toArray()
      int nUsed = usedIndices.length
      int nUsedVars = nUsed * 3
      double[] massUsed = Arrays.stream(usedIndices).
          mapToDouble({int i -> atoms[i].getAtomType().atomicWeight
          }).toArray()
      double[] xUsed = new double[nUsedVars]
      double[] x2Used = new double[nUsedVars]

      if (writeSnapshots) {
        AssemblyState origState = new AssemblyState(activeAssembly)
        forceFieldEnergy.getCoordinates(x2)
        copyCoordinates(nUsed, usedIndices, x2, x2Used)
        double[] translate = calculateTranslation(x2Used, massUsed)
        applyTranslation(x2, translate)
        forceFieldEnergy.setCoordinates(x2)
        outputFilter.writeFile(outFile, true)
        origState.revertState()
      }

      // Check which molecular assemblies to do RMSD comparisons among.
      if (!frameComparison) {
        if (dRMSD) {
          logger.info(format("\n Coordinate RMSD\n Snapshots      Original  After Translation  After Rotation    dRMSD"))
        } else if (verbose) {
          logger.info(format("\n Coordinate RMSD\n Snapshots      Original  After Translation  After Rotation"))
        }
        if (filenames.size() != 2) {
          // The first snapshot is being used for all comparisons here; therefore, snapshot = 1.
          trajectoryRMSD(systemFilter, nUsed, usedIndices, x, x2, xUsed, x2Used, massUsed, 1)
        } else {
          // Get the coordinates from the first file that was read in (i.e. an experimental PDB).
          forceFieldEnergy.getCoordinates(x)
          // The systemFilter is from the 2nd file read in, which could have multiple models.
          trajectoryRMSD(systemFilter, nUsed, usedIndices, x, x2, xUsed, x2Used, massUsed, 0)
        }
      } else {
        // Do the all vs. all comparison.
        if (storeMatrix) {
          fillDiagonals(distMatrixSize)
        }

        // Open a second copy of the system.
        MolecularAssembly[] assemblies = [potentialFunctions.open(filenames.get(0))]
        SystemFilter systemFilter2 = potentialFunctions.getFilter()

        if (dRMSD) {
          logger.info(format("\n Coordinate RMSD\n Snapshots      Original  After Translation  After Rotation    dRMSD"))
        } else if (verbose) {
          logger.info(format("\n Coordinate RMSD\n Snapshots      Original  After Translation  After Rotation"))
        }

        // Rewind the file to the first structure.
        boolean rewindFilter = true
        while (systemFilter.readNext(rewindFilter, false)) {
          rewindFilter = false
          int snapshot1 = systemFilter.getSnapshot()
          forceFieldEnergy.getCoordinates(x)

          // Compare the coordinates in x to all coordinates in the ensemble using systemFilter2.
          trajectoryRMSD(systemFilter2, nUsed, usedIndices, x, x2, xUsed, x2Used, massUsed, snapshot1)
        }
      }
    }
    return this
  }

  /**
   * Copy coordinates from the entire system to the used subset.
   *
   * @param nUsed Number of atoms used.
   * @param usedIndices Mapping from the xUsed array to its source in x.
   * @param x All atomic coordinates.
   * @param xUsed The used subset of coordinates.
   */
  private static void copyCoordinates(int nUsed, int[] usedIndices, double[] x, double[] xUsed) {
    for (int i = 0; i < nUsed; i++) {
      int index3 = 3 * usedIndices[i]
      int i3 = 3 * i
      for (int j = 0; j < 3; j++) {
        xUsed[i3 + j] = x[index3 + j]
      }
    }
  }

  /**
   * This method calculates the all versus all RMSD of a multiple model pdb/arc file.
   *
   * @param systemFilter The filter on the multiple model file.
   * @param xUsed A double array containing the xyz coordinates for multiple atoms.
   * @param x2Used A double array containing the xyz coordinates for multiple atoms.
   * @param nUsed The number of atoms that dRMSD is calculated on.
   * @param usedIndices Mapping from the xUsed array to its source in x.
   * @param x All atomic coordinates.
   * @param x2 All atomic coordinates for the second assembly.
   * @param massUsed Masses of the atoms represented in x.
   * @param snapshot1 The number of the first model being compared in the all vs. all RMSD.
   */
  void trajectoryRMSD(SystemFilter systemFilter, int nUsed, int[] usedIndices, double[] x,
      double[] x2,
      double[] xUsed, double[] x2Used, double[] massUsed, int snapshot1) {

    boolean resetFilter = true
    while (systemFilter.readNext(resetFilter, false)) {
      resetFilter = false
      int snapshot2 = systemFilter.getSnapshot()
      // Only calculate RMSD for snapshots if they aren't the same snapshot.
      // Also avoid double calculating snapshots in the matrix by only calculating the upper triangle.
      if (snapshot1 != snapshot2 && snapshot1 < snapshot2) {
        MolecularAssembly molecularAssembly2 = systemFilter.getActiveMolecularSystem()
        AssemblyState origStateB = new AssemblyState(molecularAssembly2)
        ForceFieldEnergy forceFieldEnergy2 = molecularAssembly2.getPotentialEnergy()
        forceFieldEnergy2.getCoordinates(x2)

        copyCoordinates(nUsed, usedIndices, x, xUsed)
        copyCoordinates(nUsed, usedIndices, x2, x2Used)

        double origRMSD = rmsd(xUsed, x2Used, massUsed)

        // Calculate the translation on only the used subset, but apply it to the entire structure.
        applyTranslation(x, calculateTranslation(xUsed, massUsed))
        applyTranslation(x2, calculateTranslation(x2Used, massUsed))
        // Copy the applied translation to xUsed and x2Used.
        copyCoordinates(nUsed, usedIndices, x, xUsed)
        copyCoordinates(nUsed, usedIndices, x2, x2Used)
        double translatedRMSD = rmsd(xUsed, x2Used, massUsed)

        // Calculate the rotation on only the used subset, but apply it to the entire structure.
        applyRotation(x2, calculateRotation(xUsed, x2Used, massUsed))
        // Copy the applied rotation to x2Used.
        copyCoordinates(nUsed, usedIndices, x2, x2Used)
        double rotatedRMSD = rmsd(xUsed, x2Used, massUsed)

        if (dRMSD) {
          double disRMSD = calcDRMSD(xUsed, x2Used, nUsed * 3)
          logger.info(format(
              " %6d  %6d  %7.3f            %7.3f         %7.3f  %7.3f", snapshot1, snapshot2, origRMSD,
              translatedRMSD, rotatedRMSD, disRMSD))
        } else if (verbose) {
          logger.info(format(
              " %6d  %6d  %7.3f            %7.3f         %7.3f", snapshot1, snapshot2, origRMSD, translatedRMSD,
              rotatedRMSD))
        }

        if (storeMatrix) {
          int snapshot1Index = snapshot1 - 1
          int snapshot2Index = snapshot2 - 1
          distMatrix[snapshot1Index][snapshot2Index] = rotatedRMSD
          distMatrix[snapshot2Index][snapshot1Index] = rotatedRMSD
        }

        if (writeSnapshots) {
          forceFieldEnergy.setCoordinates(x2)
          outputFilter.writeFile(outFile, true)
          origStateB.revertState()
        }
      }
    }
  }

  /**
   * Calculates the dRMSD between to sets of coordinates.
   * @param xUsed A double array containing the xyz coordinates for multiple atoms.
   * @param x2Used A double array containing the xyz coordinates for multiple atoms.
   * @param nUsed The number of atoms that dRMSD is calculated on.
   * @return A double containing the dRMSD value.
   */
  static double calcDRMSD(double[] xUsed, double[] x2Used, int nUsed) {
    double disRMSD = 0.0
    int counter = 0
    for (int i = 0; i < nUsed; i = i + 3) {
      Double3 xi = new Double3(xUsed[i], xUsed[i + 1], xUsed[i + 2])
      Double3 x2i = new Double3(x2Used[i], x2Used[i + 1], x2Used[i + 2])
      for (int j = i + 3; j < nUsed; j = j + 3) {
        Double3 xj = new Double3(xUsed[j], xUsed[j + 1], xUsed[j + 2])
        Double3 x2j = new Double3(x2Used[j], x2Used[j + 1], x2Used[j + 2])
        double dis1 = xi.sub(xj).length()
        double dis2 = x2i.sub(x2j).length()
        double diff = dis1 - dis2
        disRMSD += diff * diff
        counter++
      }
    }
    disRMSD = disRMSD / counter
    return sqrt(disRMSD)
  }

  void fillDiagonals(int size) {
    for (int i = 0; i < size; i++) {
      distMatrix[i][i] = 0.0
    }
  }

  double[][] getDistanceMatrix() {
    return distMatrix
  }
}
