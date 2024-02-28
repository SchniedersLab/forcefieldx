//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
package ffx.algorithms.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Residue
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.parsers.SystemFilter
import org.apache.commons.lang3.StringUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The Superpose script superposes molecules in an arc/multiple model pdb file (all versus all or one versus all) or in two pdb/xyz files.
 * TODO: Create a Superpose Unit Test.
 *
 * <br>
 * Usage:
 * <br>
 * ffxc Superpose [options] &lt;filename&gt;
 */
@Command(description = " Superpose frames one or two trajectory files to calculate RMSD.", name = "Superpose")
class Superpose extends AlgorithmsScript {

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  /**
   * --bb or --backboneSelection Restrict selection to backbone atoms [0 == C_ALPHA or 1 == BACKBONE]. Note C_ALPHA uses N1 or N9 for NA.
   */
  @Option(names = ['--bb', '--backboneSelection'], paramLabel = "-1", defaultValue = "-1",
      description = 'Restrict selection to backbone atoms [0 == C_ALPHA or 1 == BACKBONE]. Note C_ALPHA uses N1 or N9 for NA.')
  private int backboneSelection

  /**
   * --ih or --includeHydrogen Include hydrogen atoms.
   */
  @Option(names = ['--ih', '--includeHydrogen'], paramLabel = "false", defaultValue = "false",
      description = 'Include hydrogen atoms.')
  private static boolean includeHydrogen

  /**
   * --ss or --secondaryStructure Restrict selection using a secondary structure string.
   */
  @Option(names = ['--ss', '--secondaryStructure'], paramLabel = "", defaultValue = "",
      description = 'Restrict selection using a secondary structure string.')
  private String secondaryStructure

  /**
   * --dRMSD Calculate dRMSD in addition to RMSD.
   */
  @Option(names = ['--dRMSD'], paramLabel = "false", defaultValue = "false",
      description = 'Calculate dRMSD in addition to RMSD.')
  private boolean dRMSD

  /**
   * --ps or --printSymOp Print optimal SymOp to align structure 2 to structure 1.
   */
  @Option(names = ['--ps', '--printSymOp'], paramLabel = "false", defaultValue = "false",
          description = 'Print optimal SymOp to align input structures.')
  private static boolean printSym

  /**
   * -w or --write Write out the RMSD matrix.
   */
  @Option(names = ['-w', '--write'], paramLabel = "false", defaultValue = "false",
      description = 'Write out the RMSD matrix.')
  private static boolean write

  /**
   * -r or --restart Attempt to restart from a previously written RMSD matrix.
   */
  @Option(names = ['-r', '--restart'], paramLabel = "false", defaultValue = "false",
      description = 'Attempt to restart from a previously written RMSD matrix.')
  private static boolean restart

  /**
   * --save or --saveSnapshots Save out superposed snapshots.
   */
  @Option(names = ['--save', '--saveSnapshots'], paramLabel = "false", defaultValue = "false",
      description = 'Save superposed snapshots (only for single process jobs).')
  private boolean saveSnapshots

  /**
   * -v or --verbose Print out RMSD information.
   */
  @Option(names = ['-v', '--verbose'], paramLabel = "true", defaultValue = "true",
      description = 'Write out RMSD information.')
  private boolean verbose

  /**
   * The arguments are one or two atomic coordinate files in PDB or XYZ format.
   */
  @Parameters(arity = "1..2", paramLabel = "files",
      description = 'Atomic coordinate file(s) to compare in PDB or XYZ format.')
  List<String> filenames = null

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

    // Init the context and bind variables.
    if (!init() || filenames == null || filenames.size() < 1) {
      return this
    }

    // Turn off non-bonded terms for efficiency.
    System.setProperty("vdwterm", "false")

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filenames.get(0))
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // SystemFilter containing structures stored in file 0.
    SystemFilter baseFilter = algorithmFunctions.getFilter()

    // SystemFilter containing structures stored in file 1 (or file 0 if file 1 does not exist).
    SystemFilter targetFilter
    // Number of files to read in.
    int numFiles = filenames.size()
    boolean isSymmetric = false
    if (numFiles == 1) {
      logger.info(
          "\n Superpose will be applied between all pairs of conformations within the supplied file.\n")
      isSymmetric = true
      // If only one file is supplied, compare all structures in that file to each other.
      algorithmFunctions.openAll(filenames.get(0))
      targetFilter = algorithmFunctions.getFilter()
    } else {
      // Otherwise, compare structures from first file those in the second.
      logger.info(
          "\n Superpose will compare all conformations in the first file to all those in the second file.\n")
      algorithmFunctions.openAll(filenames.get(1))
      targetFilter = algorithmFunctions.getFilter()
    }

    // Apply active and inactive atom selections
    atomSelectionOptions.setActiveAtoms(activeAssembly)

    // Load SS definition.
    List<List<Integer>> helices = null
    List<List<Integer>> sheets = null
    boolean applySS = false
    if (!secondaryStructure.isEmpty()) {
      secondaryStructure = validateSecondaryStructurePrediction(activeAssembly)
      checkForAppropriateResidueIdentities(activeAssembly)
      String helixChar = "H"
      String sheetChar = "E"
      helices = extractSecondaryElement(secondaryStructure, helixChar, 2)
      sheets = extractSecondaryElement(secondaryStructure, sheetChar, 2)
      if (!helices.isEmpty() || !sheets.isEmpty()) {
        applySS = true
      }
    }

    // Loop over atoms and apply remaining selection flags.
    Atom[] atoms = activeAssembly.getAtomArray()
    List<Integer> selectedList = new ArrayList<>()
    for (Atom atom : atoms) {
      // Restrict active atoms to heavy atoms.
      if (!includeHydrogen && atom.isActive() && atom.isHydrogen()) {
        atom.setActive(false)
      }

      // Restrict active atoms to the SS selection.
      if (applySS && atom.isActive()) {
        if (!isSS(helices, sheets, atom)) {
          atom.setActive(false)
        }
      }

      // Restrict to the backbone atom selection.
      if (backboneSelection > -1 && atom.isActive()) {
        if (backboneSelection == 0) {
          // CALPHA (uses N1 or N9 for NA).
          if (!(atom.getName() == "CA" || atom.getName() == "N1" || atom.getName() == "N9")) {
            atom.setActive(false)
          }
        } else if (backboneSelection == 1) {
          // BACKBONE (not supported for NA)
          if (!(atom.getName() == "CA" || atom.getName() == "N" || atom.getName() == "C")) {
            atom.setActive(false)
          }
        }
      }

      // Keep only active atoms.
      if (atom.isActive()) {
        // Note that the atom index list must run from 0 to N-1.
        selectedList.add(atom.getIndex() - 1)
      }
    }

    // Indices of atoms used in alignment and RMSD calculations.
    int nSelected = selectedList.size()
    int[] usedIndices = new int[nSelected]
    for (int i = 0; i < nSelected; i++) {
      usedIndices[i] = selectedList.get(i)
    }

    // Reset all atoms to be active because the Superpose class will apply our "usedIndices" array.
    for (Atom atom : atoms) {
      atom.setActive(true)
    }

    // Compute the RMSD values.
    ffx.potential.utils.Superpose superpose =
        new ffx.potential.utils.Superpose(baseFilter, targetFilter, isSymmetric)

    // Do the superpositions.
    superpose.calculateRMSDs(usedIndices, dRMSD, verbose, restart, write, saveSnapshots, printSym)

    return this
  }

  private static boolean isSS(List<List<Integer>> helices, List<List<Integer>> sheets, Atom atom) {
    int resNum = atom.getResidueNumber() - 1
    boolean isHelix = false
    for (List<Integer> helix : helices) {
      if (resNum >= helix.get(0) && resNum <= helix.get(1)) {
        isHelix = true
        break
      }
    }
    boolean isSheet = false
    for (List<Integer> sheet : sheets) {
      if (resNum >= sheet.get(0) && resNum <= sheet.get(1)) {
        isSheet = true
        break
      }
    }
    return isHelix || isSheet
  }

  /**
   * This method determines the starting and ending indices for secondary elements of the requested type based
   * on the user-supplied secondary structure restraint predictions. The Dill Group requires that secondary
   * elements have at least three consecutive residues to be considered a secondary element.
   * @param ss A string of the secondary structure prediction.
   * @param elementType Character indicating type of secondary element being searched (helix, coil, sheet).
   * @param minNumResidues Integer minimum of consecutive secondary structure predictions
   * to create a secondary element.
   * @return ArrayList<ArrayList<Integer> > Contains starting and ending residues for each secondary element.
   */
  static List<List<Integer>> extractSecondaryElement(String ss, String elementType,
      int minNumResidues) {
    // Will hold starting and ending indices for all found secondary elements of the requested type.
    List<List<Integer>> allElements = new ArrayList<>()
    // Track of the most recent index to have a character matching the requested elementType.
    int lastMatch = 0
    // Iterates through each index in the secondary structure string.
    int i = 0
    while (i < ss.length()) {
      if (ss[i] == elementType) {
        int elementStartIndex = i
        // Set the starting index for the secondary element as soon as the value at the ith index matches the
        //  requested element type.
        for (int j = i + 1; j <= ss.length(); j++) {
          // Use the jth index to iterate through the secondary structure prediction until the end of the
          // secondary element is found.
          if (j < ss.length()) {
            if (ss[j] != elementType) {
              if (j == lastMatch + 1) {
                // If the most recent lastMatch is only one index away, then check and store the
                // starting and ending indices of the secondary element.
                i = j
                // Set i=j so that i begins searching for the next element at the end of the most recent
                // secondary element.
                int elementLength = j - elementStartIndex
                if (elementLength > minNumResidues) {
                  // If secondary element is above minimum length, store starting and ending indices
                  // of secondary element.
                  List<Integer> currentElement = new ArrayList<>()
                  currentElement.add((Integer) elementStartIndex)
                  currentElement.add((Integer) lastMatch)
                  allElements.add(currentElement)
                }
                j = ss.length() + 1
                // Since end of current secondary element has been found, exit inner j loop.
              }
            } else {
              lastMatch = j
              i++
              // If the jth index in the secondary structure string matches the requested element,
              // increment j until the end of the secondary element is found.
            }
          }
          if (j == ss.length()) {
            // Handle the case when a secondary element is at the very end of the secondary structure string.
            i = ss.length() + 1
            if (j == lastMatch + 1) {
              int elementLength = j - elementStartIndex
              if (elementLength > minNumResidues) {
                List<Integer> currentElement = new ArrayList<>()
                currentElement.add((Integer) elementStartIndex)
                currentElement.add((Integer) lastMatch)
                allElements.add(currentElement)
              }
              j = ss.length() + 1
              // Since end of current secondary element has been found, exit inner j loop.
            }
          }
        }
      } else {
        i++
        // Increment i until end of secondary structure prediction or until requested secondary element is found.
      }
    }
    return allElements
  }

  /**
   * This method validates that the user-supplied secondary structure predictions are the correct
   * length and contain the correct characters.
   *
   * @param molecularAssembly The molecular assembly.
   * @return String containing the validated and edited secondary structure restraints.
   */
  String validateSecondaryStructurePrediction(MolecularAssembly molecularAssembly) {
    // The only characters that should be present in secondary structure restraint string are 'H'
    // for helix, 'E'
    // for beta sheet and '.' for coil.
    if (!secondaryStructure.matches("^[HE.]+") && !secondaryStructure.isEmpty()) {
      logger.severe(" Secondary structure restraints may only contain characters 'H', 'E' and '.'")
    }

    int numResidues = molecularAssembly.getResidueList().size()
    int numSecondaryStructure = secondaryStructure.length()

    // Only one secondary structure restraint should exist per residue.
    if (numSecondaryStructure == 0) {
      logger.warning(
          " No secondary structure restraints have been provided. Simulation will proceed "
              + "with all residues having random coil secondary structure restraints.")
      String randomCoil = StringUtils.leftPad("", numResidues, ".")
      return randomCoil
    } else if (numSecondaryStructure < numResidues) {
      logger.warning(
          " Too few secondary structure restraints exist for number of residues present. "
              +
              "Random coil will be added to end residues without provided secondary structure restraints.")
      String extraCoil = StringUtils.rightPad(secondaryStructure, numResidues, '.')
      return extraCoil
    } else if (numSecondaryStructure == numResidues) {
      logger.info(" Secondary structure restraints will be added for all residues.")
      return secondaryStructure
    } else if (numSecondaryStructure > numResidues) {
      logger.warning(
          " Too many secondary structure restraints exist for number of residues present."
              + " Provided secondary structure restraints will be truncated.")
      String truncated = secondaryStructure.substring(0, numResidues)
      return truncated
    } else {
      logger.severe(" Secondary structure restraints or residues do not exist.")
      return null
    }
  }

  /**
   * This method checks that secondary structure assignments are appropriate for the residue
   * identity. ACE and NME residues do not have alpha carbons, so they are not compatible with the
   * alpha helix or beta sheet MELD restraints.
   *
   * @param molecularAssembly The molecular assembly.
   */
  void checkForAppropriateResidueIdentities(MolecularAssembly molecularAssembly) {
    ArrayList<Residue> residues = (ArrayList<Residue>) molecularAssembly.getResidueList()
    for (int i = 0; i < secondaryStructure.length(); i++) {
      Residue residue = residues.get(i)
      AminoAcid3 aminoAcid3 = residue.getAminoAcid3()

      String aminoAcidString = aminoAcid3.toString()
      String NMEString = AminoAcid3.NME.toString()
      String ACEString = AminoAcid3.ACE.toString()

      if (aminoAcidString == NMEString || aminoAcidString == (ACEString)) {
        char character = secondaryStructure.charAt(i)
        char H = 'H'
        char E = 'E'
        if (character == H) {
          logger.info(
              " Secondary structure was modified to accommodate non-standard amino acid residue.")
          secondaryStructure =
              secondaryStructure.substring(0, i) + '.' + secondaryStructure.substring(i + 1)
        } else if (character == E) {
          logger.info(
              " Secondary structure was modified to accommodate non-standard amino acid residue.")
          secondaryStructure =
              secondaryStructure.substring(0, i) + '.' + secondaryStructure.substring(i + 1)
        }
      }
    }
  }

}
