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
import ffx.numerics.math.RunningStatistics
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Bond
import ffx.potential.parsers.SystemFilter
import ffx.potential.utils.ProgressiveAlignmentOfCrystals
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.logging.Level

import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.concat
import static org.apache.commons.io.FilenameUtils.getBaseName
import static org.apache.commons.io.FilenameUtils.getFullPath

/**
 * Quantifies the similarity of input crystals based on progressive alignments.
 * This script is based off of the PACCOM code created by Okimasa Okada.
 *
 * @author Okimasa OKADA
 * created by Okimasa OKADA 2017/3/31
 * revised by Okimasa OKADA 2019/2/25
 * @author Aaron J. Nessler and Michael J. Schnieders
 * ported to FFX by Aaron Nessler and Micheal Schnieders 2020
 * revised by Aaron Nessler and Michael Schnieders 2021
 * <br>
 * Usage:
 * <br>
 * ffxc SuperposeCrystals &lt;filename&gt; &lt;filename&gt;
 */
@Command(description = " Determine the RMSD for crystal polymorphs using the Progressive Alignment of Crystals (PAC) algorithm.",
    name = "SuperposeCrystals")
class SuperposeCrystals extends AlgorithmsScript {

  /**
   * --ac or --alchemicalAtoms sets atoms unique for both crystals, as comma-separated hyphenated
   * ranges or singletons.
   */
  @Option(names = ["--ac", "--alchemicalAtoms"], paramLabel = "", defaultValue = "",
          description = "Atom indices to be excluded from both crystals (e.g. 1-24,32-65). Use if molecular identity and atom ordering are the same for both crystals (otherwise use \"--ac1\" and \"--ac2\".")
  private String excludeAtoms = ""

  /**
   * --ac1 or --alchemicalAtoms1 sets atoms unique to the base crystal, as comma-separated hyphenated
   * ranges or singletons.
   */
  @Option(names = ["--ac1", "--alchemicalAtoms1"], paramLabel = "", defaultValue = "",
          description = "Atom indices to be excluded in the first crystal (e.g. 1-24,32-65).")
  private String excludeAtomsA = ""

  /**
   * --ac2 or --alchemicalAtoms2 sets atoms unique to the target crystal, as comma-separated hyphenated
   * ranges or singletons.
   */
  @Option(names = ["--ac2", "--alchemicalAtoms2"], paramLabel = "", defaultValue = "",
          description = "Atom indices to be excluded in the second crystal (e.g. 1-24,32-65).")
  private String excludeAtomsB = ""

  /**
   * --ca or --carbonAlphas Consider only alpha carbons for proteins.
   */
  @Option(names = ['--ca', '--carbonAlphas'], paramLabel = "false", defaultValue = "false",
          description = 'Consider only alpha carbons for proteins.')
  private static boolean alphaCarbons

  /**
   * --cd or --createDirectories Create subdirectories for FE simulations.
   */
  @Option(names = ['--cd', '--createDirectories'], paramLabel = "false", defaultValue = "false",
          description = 'Create subdirectories for free energy simulations.')
  private static boolean createFE

  /**
   * --as or --autoSymmetry Automatic generation of symmetry operators.
   */
  @Option(names = ['--as', '--autoSymmetry'], paramLabel = "false", defaultValue = "false",
          description = 'Automatically generate symmetry operators.')
  private static boolean autoSym

  /**
   * --st or --symTolerance Tolerance used to add small atom groups for automatic symmetry matching.
   */
  @Option(names = ['--st','--symTolerance'], paramLabel = "0.5", defaultValue = "0.5",
          description = 'Add atoms to torsion based symmetries that are beneath this cutoff.')
  private static double symTolerance

  /**
   * --ws or --writeSym Write symmetry operator to Properties/Key file.
   */
  @Option(names = ['--ws', '--writeSym'], paramLabel = "false", defaultValue = "false",
          description = 'Write the created symmetry operators in the corresponding Properties/Key file.')
  private static boolean writeSym

  /**
   * --gc or --gyrationComponents Display components for radius of gyration for final clusters.
   */
  @Option(names = ['--gc', '--gyrationComponents'], paramLabel = "false", defaultValue = "false",
          description = 'Display components for radius of gyration for final clusters.')
  private static boolean gyrationComponents

  /**
   * --ht or --hitTolerance Tolerance to determine if a comparison should be counted as a "hit".
   */
  @Option(names = ['--ht', '--hitTolerance'], paramLabel = '-1.0', defaultValue = '-1.0',
          description = "Sum comparisons that attain a value lower than this tolerance.")
  private double hitTol

  /**
   * --if or --inflationFactor Inflation factor used to determine replicates expansion.
   */
  @Option(names = ['--if', '--inflationFactor'], paramLabel = '5.0', defaultValue = '5.0',
          description = 'Inflation factor used to determine replicates expansion (IF * nAU in replicates).')
  private double inflationFactor

  /**
   * --ih or --includeHydrogen Include hydrogen atoms.
   */
  @Option(names = ['--ih', '--includeHydrogen'], paramLabel = "false", defaultValue = "false",
          description = 'Include hydrogen atoms.')
  private boolean includeHydrogen

  /**
   * --in or --inertia Display moments of inertia for final clusters.
   */
  @Option(names = ['--in', '--inertia'], paramLabel = "false", defaultValue = "false",
          description = 'Display moments of inertia for final clusters.')
  private static boolean inertia

  /**
   * -l or --linkage Single (0), Average (1), or Complete (2) coordinate linkage for molecule prioritization.
   */
  @Option(names = ['-l', '--linkage'], paramLabel = '1', defaultValue = '1',
          description = 'Single (0), Average (1), or Complete (2) coordinate linkage for molecule prioritization.')
  private int linkage

  /**
   * --lm or --lowMemory Slower comparisons, but reduces memory usage.
   */
  @Option(names = ['--lm', '--lowMemory'], paramLabel = "false", defaultValue = "false",
          description = 'Reduce memory usage at the cost of efficiency.')
  private static boolean lowMemory

  /**
   * --mt or --matchTolerance Tolerance to determine if two AUs are different.
   */
  @Option(names = ['--mt', '--moleculeTolerance'], paramLabel = '0.0015', defaultValue = '0.0015',
          description = "Tolerance to determine if two AUs are different.")
  private double matchTol

  /**
   * --mw or --massWeighted Use mass-weighted atomic coordinates for alignment.
   */
  @Option(names = ['--mw', '--massWeighted'], paramLabel = "false", defaultValue = "false",
          description = 'Use mass-weighted atomic coordinates for alignment.')
  private static boolean massWeighted

  /**
   * --na or --numAU AUs in the RMSD.
   */
  @Option(names = ['--na', '--numAU'], paramLabel = '20', defaultValue = '20',
          description = 'Number of asymmetric units included to calculate the RMSD.')
  private int numAU

  /**
   * --pc or --prioritizeCrystals Prioritize the crystals being compared based on high density (0), low density (1), or file order (2).
   */
  @Option(names = ['--pc', '--prioritizeCrystals'], paramLabel = '0', defaultValue = '0',
          description = 'Prioritize crystals based on high density (0), low density (1), or file order (2).')
  private int crystalPriority

  /**
   * --ps or --printSymOp Print optimal SymOp to align crystal 2 to crystal 1 (atoms deviating above this value are highlighted).
   */
  @Option(names = ['--ps', '--printSymOp'], paramLabel = "-1.0", defaultValue = "-1.0",
          description = 'Print optimal SymOp to align input crystals (print out atom deviations above value).')
  private static double printSym

  /**
   * -r or --restart Restart from a previously written RMSD matrix (if one exists).
   */
  @Option(names = ['-r', '--restart'], paramLabel = "false", defaultValue = "false",
          description = 'Restart from a previously written RMSD matrix (if one exists).')
  private static boolean restart

  /**
   * -s or --save Save files for the superposed crystals that align below tolerance.
   */
  @Option(names = ['-s','--save'], paramLabel = "-1.0", defaultValue = "-1.0",
          description = 'Save structures less than or equal to this cutoff.')
  private static double save

  /**
   * --sc --saveClusters Save files for the superposed crystal clusters.
   */
  @Option(names = ['--sc', '--saveClusters'], paramLabel = "0", defaultValue = "0",
          description = 'Save files for the superposed crystal clusters (1=PDB, 2=XYZ).')
  private static int saveClusters

  /**
   * --sm or --saveMachineLearning Save out PDB and CSV for machine learning.
   */
  @Option(names = ['--sm', '--saveMachineLearning'], paramLabel = "false", defaultValue = "false",
          description = 'Final structures for each comparison will be written out with RMSD in a CSV.')
  private static boolean machineLearning

  /**
   * --th or --thorough Compare all unique AUs between each crystal.
   */
  @Option(names = ['--th', '--thorough'], paramLabel = "false", defaultValue = "false",
          description = 'More intensive, less efficient version of PAC.')
  private static boolean thorough

  /**
   * -w or --write Write out the PAC RMSD matrix.
   */
  @Option(names = ['-w', '--write'], paramLabel = "false", defaultValue = "false",
          description = 'Write out the PAC RMSD matrix.')
  private static boolean write

  /**
   * --zp or --zPrime Z' for both crystals (same).
   */
  @Option(names = ['--zp', '--zPrime'], paramLabel = '-1', defaultValue = '-1',
          description = "Z'' for both crystals (assumes same value).")
  private int zPrime

  /**
   * --zp1 or --zPrime1 Z' for crystal 1 (default to autodetect).
   */
  @Option(names = ['--zp1', '--zPrime1'], paramLabel = '-1', defaultValue = '-1',
          description = "Z'' for crystal 1 (default will try to autodetect).")
  private int zPrime1

  /**
   * --zp2 or --zPrime2 Z' for crystal 2 (default to autodetect).
   */
  @Option(names = ['--zp2', '--zPrime2'], paramLabel = '-1', defaultValue = '-1',
          description = "Z'' for crystal 2 (default will try to autodetect).")
  private int zPrime2

  /**
   * The final argument(s) should be two or more filenames (same file twice if comparing same structures).
   */
  @Parameters(arity = "1..2", paramLabel = "files",
          description = 'Atomic coordinate file(s) to compare in XYZ format.')
  List<String> filenames = null

  /**
   * List of atoms that have already been included for automatic symmetry generation.
   */
  private List<Atom> previouslyIncluded = new ArrayList<>();

  /**
   * List of indicies that may be added to automatic symmetry generation.
   */
  private List<Integer> queueIndices = new ArrayList<>();

  /**
   * CrystalSuperpose Test requires a public variable containing observables to test.
   */
  public RunningStatistics runningStatistics

  /**
   * Symmetry operators applied to
   */
  public final StringBuilder symOpsA = new StringBuilder("");
  /**
   * Symmetry operators applied to
   */
  public final StringBuilder symOpsB = new StringBuilder("");

  /**
   * CrystalSuperpose Constructor.
   */
  SuperposeCrystals() {
    this(new Binding())
  }

  /**
   * CrystalSuperpose Constructor.
   * @param binding Groovy Binding to use.
   */
  SuperposeCrystals(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  SuperposeCrystals run() {

    // Turn off non-bonded interactions for speed.
    System.setProperty("vdwterm", "false")

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Ensure file exists.
    if (filenames == null) {
      logger.info(helpString())
      return this
    }

    // SystemFilter containing structures stored in file 0.
    algorithmFunctions.openAll(filenames.get(0))
    SystemFilter baseFilter = algorithmFunctions.getFilter()
    // SystemFilter containing structures stored in file 1 (or file 0 if file 1 does not exist).
    SystemFilter targetFilter

    // Number of files to read in.
    boolean isSymmetric = false
    int numFiles = filenames.size()
    if (numFiles == 1) {
      logger.info(
              "\n PAC will be applied between all pairs of conformations within the supplied file.\n")
      isSymmetric = true
      // If only one file is supplied, compare all structures in that file to each other.
      algorithmFunctions.openAll(filenames.get(0))
      targetFilter = algorithmFunctions.getFilter()
    } else {
      // Otherwise, compare structures from first file those in the second.
      logger.info(
              "\n PAC will compare all conformations in the first file to all those in the second file.\n")
      algorithmFunctions.openAll(filenames.get(1))
      targetFilter = algorithmFunctions.getFilter()
      String filenameA = getBaseName(filenames.get(0));
      String filenameB = getBaseName(filenames.get(1));
      symOpsA.append(format("\n Add the following symop to properties/key file of %s (Apply on %s to generate %s):\nsymop", filenameA, filenameB, filenameA));
      symOpsB.append(format("\n\n Add the following symop to properties/key file of %s (Apply on %s to generate %s):\nsymop", filenameB, filenameA, filenameB));
    }

    // Define the filename to use for the RMSD values.
    String filename = filenames.get(0)
    String pacFilename = concat(getFullPath(filename), getBaseName(filename) + ".dst")

    if (zPrime > 0 && zPrime % 1 == 0) {
      zPrime1 = zPrime;
      zPrime2 = zPrime;
    }

    // Write sym and automatic sym generation both depend on calculating symmetry operators, set --ps if not already.
    if (printSym < 0 && writeSym || autoSym ) {
      logger.info("\n Printing distance for atoms greater than large distance of 10.0 Å");
      printSym = 10.0;
    }

    if(createFE){
      // TODO see about caching entire assemblies rather than pieces.
      // Need assembly to generate subdirectories. Caching ignores assemblies.
      lowMemory = true;
    }

    // The autoSym flag is heuristic based... May need updated as more structures are tested.
    //  AutoSym loops through and groups atoms into "small" groups 3+ atoms (important for multipole orientation).
    //  Next loops through small groups and creates final groups according to symTolerance cutoff.
    if(autoSym) {
      // Reduce logging to final groups if using default logging (small group PAC produces a lot of logging).
      boolean loggerLevelChanged = false
      if(logger.isLoggable(Level.INFO) && !logger.isLoggable(Level.FINE)){
        loggerLevelChanged = true;
        logger.setLevel(Level.WARNING);
      }
      // List to contain atoms that have already been grouped. Therefore, skip for future groups.
      final Atom[] atoms = baseFilter.getActiveMolecularSystem().getAtomArray();
      final int nAtoms = atoms.length;
      // Need 3 atoms to define multipole axes of sym groups, if less, combine groups until sufficient.
      // Retain list of small groups to combine into larger groups later.
      List<List<Atom>> smallAtomGroups = new ArrayList<ArrayList<Atom>>();
      // Record small group RMSD to know if they can be combined (if RMSD > tolerance, then should be independent group).
      List<Double> smallGroupsRMSD = new ArrayList<Double>();
      // Loop through atoms in system.
      for (Atom a : atoms) {
        // Skip if already included.
        if (previouslyIncluded.contains(a)) {
          continue;
        }
        // Clear queued indices if there are already any atoms included.
        queueIndices.clear();
        // Add this atom to list (recursive method).
        addAtomIndex(a);

        ArrayList<Atom> queueAtoms = new ArrayList<Atom>();
        for (Integer index : queueIndices) {
          // Indices are human readable ergo (-1).
          queueAtoms.add(atoms[index - 1]);
        }

        // This is a kludgy fail safe... Essentially throw inappropriately grouped atoms with a bonded atom's group.
        if (queueAtoms.size() < 3) {
          boolean groupFound = false;
          // Loop through group that is too small to satisfy multipoles.
          for (Atom atom : queueAtoms) {
            Bond[] bonds = atom.getBonds()
            // Find an atom bonded to atom in this group.
            for (Bond bond : bonds) {
              Atom atomGrouped = bond.get1_2(atom);
              // Find small group corresponding to bonded atom.
              for (List<Atom> atomList : smallAtomGroups) {
                if (atomList.contains(atomGrouped)) {
                  groupFound = true;
                  // Add all atoms from group that is too small to bonded atom's group.
                  for (Atom atomNeedingGroup : queueAtoms) {
                    logger.warning(" Atom needs group: " + atomNeedingGroup.toString());
                    atomList.add(atomNeedingGroup);
                  }
                }
                if (groupFound) {
                  break;
                }
              }
              if (groupFound) {
                break;
              }
            }
            if (groupFound) {
              break;
            }
          }
        } else {
          // This is the normal/expected behavior assuming correct group size.
          smallAtomGroups.add(queueAtoms);
        }
      }

      // Check for atoms that have not been included (this code should not be needed... Seems like -N- or -O- ).
      for (int i = 0; i < nAtoms; i++) {
        Atom atom = atoms[i]
        boolean found = false;
        for(List<Atom> list : smallAtomGroups) {
          if(list.contains(atom)){
            found = true;
            break;
          }
        }
        if(!found){
          logger.warning(" ATOM NOT INCLUDED: " + atom.toString());
          boolean groupFound = false;
          Bond[] bonds = atom.getBonds()
          // Find an atom bonded to atom.
          for (Bond bond : bonds) {
            Atom atomGrouped = bond.get1_2(atom);
            // Find small group corresponding to bonded atom.
            for (List<Atom> atomList : smallAtomGroups) {
              if (atomList.contains(atomGrouped)) {
                groupFound = true;
                // Add atom to bonded atom's group.
                atomList.add(atom);
              }
              if (groupFound) {
                break;
              }
            }
            if (groupFound) {
              break;
            }
          }
        }
      }

      // Atoms have been grouped into sufficient groups of minimal size.
      StringBuilder sb = new StringBuilder(" Small Groups:\n");
      for(ArrayList<Atom> group: smallAtomGroups){
        // Run PAC to assess small group selection.
        ArrayList<Integer> includeIndices = new ArrayList<>();
        sb.append(" New Group: ");
        for(Atom atom: group){
          sb.append(atom.getIndex() + ", ");
          includeIndices.add(atom.getIndex());
        }
        sb.append("\n");

        excludeAtoms = "";

        for (int i = 0; i < nAtoms; i++) {
          int ind = i + 1;
          if (!includeIndices.contains(ind) && excludeAtoms.length() == 0) {
            excludeAtoms = ind;
          } else if (!includeIndices.contains(ind)) {
            excludeAtoms += "," + ind;
          }
        }
        // Apply atom selections
        if (excludeAtoms == null && excludeAtoms.isEmpty()) {
          continue;
        } else {
          excludeAtomsA = excludeAtoms
          excludeAtomsB = excludeAtoms
        }
        if (logger.isLoggable(Level.FINE)) {
          logger.fine("A: " + excludeAtomsA);
          logger.fine("B: " + excludeAtomsB);
        }

        // Compare structures in baseFilter and targetFilter.
        ProgressiveAlignmentOfCrystals pac = new ProgressiveAlignmentOfCrystals(baseFilter, targetFilter,
                isSymmetric)
        runningStatistics =
                pac.comparisons(numAU, inflationFactor, matchTol, hitTol, zPrime1, zPrime2, excludeAtomsA, excludeAtomsB,
                        alphaCarbons, includeHydrogen, massWeighted, crystalPriority, thorough, saveClusters, save,
                        restart, write, machineLearning, inertia, gyrationComponents, linkage, printSym,
                        lowMemory, createFE, false, pacFilename, new StringBuilder(""), new StringBuilder(""))
        double min = runningStatistics.min;
        smallGroupsRMSD.add(min)
      }
      logger.info(sb.toString());
      int numGroups = smallAtomGroups.size();
      if(logger.isLoggable(Level.FINE)){
        logger.fine(format(" Number of minimal groups: %3d", numGroups));
      }
      int[][] smallIndexGroups = new int[numGroups][];
      // Small atom groups are set. Create jagged array containing indices.
      for(int i = 0; i < numGroups; i++){
        // Obtain atoms for this group.
        ArrayList<Atom> currentAtoms = smallAtomGroups.get(i);
        // Determine size of this row.
        int currentSize = currentAtoms.size()
        smallIndexGroups[i] = new int[currentSize];
        // Add indices of atoms to this row.
        for(int j = 0; j < currentSize; j++){
          smallIndexGroups[i][j] = currentAtoms.get(j).getIndex();
        }
      }
      // Object to hold the final groups to be symmetry operated.
      ArrayList<ArrayList<Integer>> finalGroups = new ArrayList<>();
      // Atoms that will be included in this group.
      ArrayList<Integer> acceptedAtoms = new ArrayList<>();
      // Atoms to attempt to include in this group.
      ArrayList<Integer> queuedIndices = new ArrayList<>();
      // Groups that have already been mapped.
      ArrayList<Integer> usedGroups = new ArrayList<>();
      // Loop through each of the small atom groups.
      for(int i = 0; i < numGroups; i++) {
        ArrayList<Atom> currentAtomGroup = smallAtomGroups.get(i);
        ArrayList<Integer> currentIndexGroup = smallIndexGroups[i];
        int numAtomsInGroup = currentAtomGroup.size();
        // Reset object lists for new group to be added to finalGroups.
        acceptedAtoms.clear();
        queuedIndices.clear();
        if(logger.isLoggable(Level.FINE)) {
          logger.fine(format(" Minimal group number: %3d Size: %3d RMSD: %9.3f", i, numAtomsInGroup, smallGroupsRMSD.get(i)))
        }
        if (!usedGroups.contains(i)) {
          // Add current group under consideration to accepted atoms.
          for (Integer index : currentIndexGroup) {
            acceptedAtoms.add(index);
          }
          usedGroups.add(i);
          // Groups already above symTolerance cannot be added to other groups and are included outright.
          if (smallGroupsRMSD.get(i) > symTolerance) {
            if(logger.isLoggable(Level.FINE)) {
              logger.fine(format(" Run (group) %3d added outright to final.", i));
            }
          }else{
            for (int j = 0; j < acceptedAtoms.size(); j++){
              Atom a1 = atoms[acceptedAtoms.get(j) - 1];
              if(logger.isLoggable(Level.FINE)) {
                logger.fine(format(" Current group: %3d Atom in Group: %3d of %3d (was %3d==%3d) Atom Index: %3d", i, j, acceptedAtoms.size(), numAtomsInGroup, currentIndexGroup.size(), a1.getIndex()));
              }
              Bond[] b1 = a1.getBonds();
              // Traverse bonds to find adjacent small groups.
              for (Bond b : b1) {
                Atom a2 = b.get1_2(a1);
                int a2Index = a2.getIndex();
                // If atom is encountered that is not currently in this group, evaluate to see if it should be added.
                if (!currentIndexGroup.contains(a2Index)) {
                  // Find the group that this new atom belongs.
                  for (int k = 0; k < numGroups; k++) {
                    ArrayList<Integer> currentIndexGroup2 = smallIndexGroups[k];
                    // Group cannot be added to two separate sym ops.
                    // If sym is greater than tolerance it cannot be added.
                    // If atom is not included in this group, go to next.
                    if(!usedGroups.contains(k) && smallGroupsRMSD[k] < symTolerance && currentIndexGroup2.contains(a2Index)) {
                      // Attempt to add group containing bonded atom to current group.
                      for (Integer ind : currentIndexGroup2) {
                        queuedIndices.add(ind);
                      }
                      // Create exclusion strings to isolate desired atoms in PAC.
                      excludeAtoms = "";
                      for (int l = 0; l < nAtoms; l++) {
                        int ind = l + 1;
                        // Any atoms that are not in the accepted or queued lists should be excluded.
                        boolean exclude = !acceptedAtoms.contains(ind) && !queuedIndices.contains(ind)
                        if (exclude) {
                          // Format input string to PAC.
                          if (excludeAtoms.length() == 0) {
                            excludeAtoms = ind;
                          } else {
                            excludeAtoms += "," + ind;
                          }
                        }
                      }
                      // Apply atom selections
                      excludeAtomsA = excludeAtoms
                      excludeAtomsB = excludeAtoms

                      if(logger.isLoggable(Level.FINE)) {
                        logger.fine(" Excluded atoms A: " + excludeAtomsA);
                        logger.fine(" Excluded atoms B: " + excludeAtomsB);
                      }

                      // Compare structures in baseFilter and targetFilter.
                      ProgressiveAlignmentOfCrystals pac = new ProgressiveAlignmentOfCrystals(baseFilter, targetFilter,
                              isSymmetric)
                      runningStatistics =
                              pac.comparisons(numAU, inflationFactor, matchTol, hitTol, zPrime1, zPrime2, excludeAtomsA, excludeAtomsB,
                                      alphaCarbons, includeHydrogen, massWeighted, crystalPriority, thorough, saveClusters, save,
                                      restart, write, machineLearning, inertia, gyrationComponents, linkage, printSym,
                                      lowMemory, createFE, false, pacFilename, new StringBuilder(""), new StringBuilder(""))
                      double min = runningStatistics.min;
                      if (min > symTolerance) {
                        // Group should not be added as it increased RMSD beyond tolerance.
                        if(logger.isLoggable(Level.FINE)) {
                          logger.fine(format(" Run %3d Group %3d failed with RMSD %9.3f > tolerance (%9.3f).", i, k, min, symTolerance));
                        }
                        queuedIndices.clear();
                      } else {
                        // Satisfied RMSD cutoff, therefore add to accepted atoms.
                        for (Integer value : queuedIndices) {
                          if (!acceptedAtoms.contains(value)) {
                            acceptedAtoms.add(value);
                          }
                        }
                        if(logger.isLoggable(Level.FINE)) {
                          logger.fine(format(" Run %3d Group %3d added to accepted with RMSD: %9.3f.", i, k, min));
                        }
                        usedGroups.add(k);
                      }
                    }
                  }
                }
              }
            }
          }
          finalGroups.add((ArrayList<Integer>) acceptedAtoms.clone());
          if(logger.isLoggable(Level.FINE)) {
            logger.fine(format(" Run %3d accepted atoms (size: %3d) added to finalGroups (size: %3d)", i, acceptedAtoms.size(), finalGroups.size()))
          }
        }
      }

      if(logger.isLoggable(Level.FINE)){
        logger.fine(format(" Final Num Groups: %3d", finalGroups.size()))
      }

      // Small atom groups finished. Update logging for final groups.
      if(loggerLevelChanged){
        logger.setLevel(Level.INFO);
      }
      for (ArrayList<Integer> group : finalGroups) {
        int groupSize = group.size();
        if(logger.isLoggable(Level.FINE)){
          logger.fine(format(" Final Size: %3d", groupSize))
        }
        // Skip groups with size of zero. THIS SHOULD NEVER OCCUR.
        if(groupSize <= 0){
          logger.warning(" Final group of size zero was encountered. Check selection logic.")
          continue;
        }
        queuedIndices.clear();
        for(Integer ind : group){
          queuedIndices.add(ind);
        }
        excludeAtoms = "";
        for (int l = 0; l < nAtoms; l++) {
          int ind = l + 1;
          if (!queuedIndices.contains(ind) && excludeAtoms.length() == 0) {
            excludeAtoms = ind;
          } else if (!queuedIndices.contains(ind)) {
            excludeAtoms += "," + ind;
          }
        }
        // Apply atom selections
        excludeAtomsA = excludeAtoms
        excludeAtomsB = excludeAtoms
        if (logger.isLoggable(Level.FINE)) {
          logger.fine("A: " + excludeAtomsA);
          logger.fine("B: " + excludeAtomsB);
        }

        // Compare structures in baseFilter and targetFilter.
        ProgressiveAlignmentOfCrystals pac = new ProgressiveAlignmentOfCrystals(baseFilter, targetFilter,
                isSymmetric)
        runningStatistics =
                pac.comparisons(numAU, inflationFactor, matchTol, hitTol, zPrime1, zPrime2, excludeAtomsA, excludeAtomsB,
                        alphaCarbons, includeHydrogen, massWeighted, crystalPriority, thorough, saveClusters, save,
                        restart, write, machineLearning, inertia, gyrationComponents, linkage, printSym,
                        lowMemory, createFE, writeSym, pacFilename, symOpsA, symOpsB)
      }
    } else {
      // Apply atom selections
      if (excludeAtoms != null && !excludeAtoms.isEmpty()) {
        excludeAtomsA = excludeAtoms
        excludeAtomsB = excludeAtoms
      }
      // Compare structures in baseFilter and targetFilter.
      ProgressiveAlignmentOfCrystals pac = new ProgressiveAlignmentOfCrystals(baseFilter, targetFilter,
              isSymmetric)
      runningStatistics =
              pac.comparisons(numAU, inflationFactor, matchTol, hitTol, zPrime1, zPrime2, excludeAtomsA, excludeAtomsB,
                      alphaCarbons, includeHydrogen, massWeighted, crystalPriority, thorough, saveClusters, save,
                      restart, write, machineLearning, inertia, gyrationComponents, linkage, printSym,
                      lowMemory, createFE, writeSym, pacFilename, symOpsA, symOpsB)
    }

    // Print out symmetry operators
    if (printSym >= 0.0) {
      int finalCharA = symOpsA.size() - 2
      if (symOpsA.substring(finalCharA) == '\\\n') {
        symOpsA.setLength(finalCharA)
      }
      int finalCharB = symOpsB.size() - 2
      if (symOpsB.substring(finalCharB) == '\\\n') {
        symOpsB.setLength(finalCharB)
      }
      if (symOpsA.contains("\\") && symOpsB.contains("\\")) {
        logger.info(symOpsA.toString() + symOpsB.toString());
      } else if (logger.isLoggable(Level.FINE)) {
        logger.fine(" Following generated SymOp strings did not contain symmetry operators." + symOpsA.toString() + symOpsB.toString())
      }
    }

    // Close readers.
    baseFilter.closeReader();
    targetFilter.closeReader();
    return this
  }

  /**
   * Search of bonded atoms to those added to the queue list to look for additional inclusions.
   * @param aAdded Atom to be added to the current group.
   * @return
   */
  private void addAtomIndex(Atom aAdded){
    if(!previouslyIncluded.contains(aAdded)) {
      queueIndices.add(aAdded.getIndex());
      previouslyIncluded.add(aAdded)
      int numBondsAdded = aAdded.getNumBonds();
      int numHAdded = aAdded.getNumberOfBondedHydrogen();
      int nonHbonds = numBondsAdded - numHAdded;

      Bond[] bonds4Added = aAdded.getBonds();
      for (Bond b : bonds4Added) {
        Atom aBound2Added = b.get1_2(aAdded);
        if (nonHbonds == 1) {
          // If only bond to group add (e.g., -O-CH3, -O-NH2).
          addAtomIndex(aBound2Added)
          for (Bond bAdded : bonds4Added) {
            Atom a3 = bAdded.get1_2(aAdded);
            if (a3.isHydrogen()) {
              addAtomIndex(a3)
            }
          }
        }

        int numBondsB2A = aBound2Added.getNumBonds();
        int numHB2A = aBound2Added.getNumberOfBondedHydrogen();
        int nonHB2A = numBondsB2A - numHB2A;
        Bond[] bondsB2A = aBound2Added.getBonds();
        if (nonHB2A == 1) {
          // If only bond to group add (e.g., -O-CH3, -O-NH2).
          addAtomIndex(aBound2Added);
          for (Bond b2a : bondsB2A) {
            Atom aBond2B2A = b2a.get1_2(aBound2Added);
            if (aBond2B2A.isHydrogen()) {
              addAtomIndex(aBond2B2A);
            }
          }
        } else if (numBondsAdded == 1 || numBondsB2A == 1) {
          // If second atom only has one bond and it bonds to 'a' add to list (e.g., hydrogen and halogen)
          addAtomIndex(aBound2Added)
        } else if (nonHB2A == 2) {
          // If second atom is only connected to this and one other group add (e.g., )
          addAtomIndex(aBound2Added)
          for (Bond b2a : bondsB2A) {
            Atom aBond2B2A = b2a.get1_2(aBound2Added);
            if (aBond2B2A.isHydrogen()) {
              addAtomIndex(aBond2B2A);
            }
          }
        } else {
          // TODO handle higher order bonds and aromaticity.
          // Look for 1_3 atoms that solely attach to a2 (e.g., -C=O-NH2, -C=0-CH3)
          for (Bond b2a : bondsB2A) {
            Atom aBond2B2A = b2a.get1_2(aBound2Added);
            if (!aBond2B2A.isHydrogen() && nonHbonds == 1) {
              int numBonds2B2A = aBond2B2A.getNumBonds();
              if ((numBonds2B2A == 1) && (aBond2B2A.isBonded(aBound2Added) || aBond2B2A.isBonded(aAdded))) {
                addAtomIndex(aBond2B2A);
              }
            }
          }
        }
      }
    }
  }
}