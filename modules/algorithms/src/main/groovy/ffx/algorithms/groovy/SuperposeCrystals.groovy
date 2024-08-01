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
   * --bs or --bruteSymmetry Brute force symmetry operator creation.
   */
  @Option(names = ['--bs', '--bruteSymmetry'], paramLabel = "false", defaultValue = "false",
          description = 'Brute force symmetry operator creation.')
  private static boolean bruteSym

  /**
   * --as or --appendSym Append symmetry operator.
   */
  @Option(names = ['--as', '--appendSym'], paramLabel = "false", defaultValue = "false",
          description = 'Append the created symmetry operators.')
  private static boolean appendSym

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
   * --ps or --printSymOp Print optimal SymOp to align crystal 2 to crystal 1.
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
   * --save Save files for the superposed crystals.
   */
  @Option(names = ['--save'], paramLabel = "-1.0", defaultValue = "-1.0",
          description = 'Save structures less than or equal to this cutoff.')
  private static double save

  /**
   * --sc --saveClusters Save files for the superposed crystals.
   */
  @Option(names = ['--sc', '--saveClusters'], paramLabel = "0", defaultValue = "0",
          description = 'Save files for the superposed crystals (1=PDB, 2=XYZ).')
  private static int saveClusters

  /**
   * --sm or --saveMachineLearning Save out PDB and CSV for machine learning.
   */
  @Option(names = ['--sm', '--saveMachineLearning'], paramLabel = "false", defaultValue = "false",
          description = 'Final structures for each comparison will be written out with RMSD in a CSV.')
  private static boolean machineLearning

  /**
   * --st or --strict Compare all unique AUs between each crystal.
   */
  @Option(names = ['--st', '--strict'], paramLabel = "false", defaultValue = "false",
          description = 'More intensive, less efficient version of PAC.')
  private static boolean strict

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
   * CrystalSuperpose Test requires a public variable containing observables to test.
   */
  public RunningStatistics runningStatistics

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
    }

    // Define the filename to use for the RMSD values.
    String filename = filenames.get(0)
    String pacFilename = concat(getFullPath(filename), getBaseName(filename) + ".dst")

    if (zPrime > 0 && zPrime % 1 == 0) {
      zPrime1 = zPrime;
      zPrime2 = zPrime;
    }

    if (appendSym && printSym < 0) {
      logger.info(" Printing distance for atoms greater than default of 1.0 Å");
      printSym = 1.0;
    }
    if(createFE){
      // TODO see about caching entire assemblies rather than pieces.
      // Need assembly to generate subdirectories. Caching ignores assemblies.
      lowMemory = true;
    }
    if(bruteSym) {
      List<Atom> previouslyIncluded = new ArrayList<>();
      final Atom[] atoms = baseFilter.getActiveMolecularSystem().getAtomArray();
      for (Atom a : atoms) {
        if (previouslyIncluded.contains(a)) {
          continue;
        }
        // Collect information used for decision making.
        Bond[] bondsa = a.getBonds();
        int numBondsA = a.getNumBonds();
        int numHA = a.getNumberOfBondedHydrogen();
        int nonHbonds = numBondsA - numHA;
        List<Integer> includeAtoms = new ArrayList<>();
        List<Integer> queueAtoms = new ArrayList<>();
        queueAtoms.add(a.getIndex());
        previouslyIncluded.add(a);
        for (Atom a2 : atoms) {
          // If already included, skip this atom. Compare returns 0 when equal.
          if (previouslyIncluded.contains(a2) || !a2.isBonded(a)) {
            continue;
          }
          // Collect information used for decision making.
          int numBondsA2 = a2.getNumBonds();;
          int numHA2 = a2.getNumberOfBondedHydrogen();
          int nonHbonds2 = numBondsA2 - numHA2;
          Bond[] bondsa2 = a2.getBonds();
          // Decide which atoms should be grouped together.
          if (nonHbonds2 == 1) {
            // If only bond to group add (e.g., -O-CH3, -O-NH2).
            queueAtoms.add(a2.getIndex());
            previouslyIncluded.add(a2);
            for(Bond b: bondsa2){
              Atom a3 = b.get1_2(a2);
              if(previouslyIncluded.contains(a3)){
                continue;
              }else if(a3.isHydrogen()){
                queueAtoms.add(a3.getIndex());
                previouslyIncluded.add(a3);
              }
            }
          } else if (nonHbonds == 1) {
            // If only bond to group add (e.g., -O-CH3, -O-NH2).
            queueAtoms.add(a2.getIndex());
            previouslyIncluded.add(a2);
            for (Bond b: bondsa) {
              Atom a3 = b.get1_2(a);
              if (previouslyIncluded.contains(a3)) {
                continue;
              } else if (a3.isHydrogen()) {
                queueAtoms.add(a3.getIndex());
                previouslyIncluded.add(a3);
              }
            }
          } else if (numBondsA == 1 || numBondsA2 == 1) {
            // If second atom only has one bond and it bonds to 'a' add to list (e.g., hydrogen and halogen)
            queueAtoms.add(a2.getIndex());
            previouslyIncluded.add(a2);
          } else if (nonHbonds2 == 2) {
            // If second atom is only connected to this and one other group add (e.g., )
            queueAtoms.add(a2.getIndex());
            previouslyIncluded.add(a2);
            for(Bond b: bondsa2){
              Atom a3 = b.get1_2(a2);
              if(previouslyIncluded.contains(a3)){
                continue;
              }else if(a3.isHydrogen()){
                queueAtoms.add(a3.getIndex());
                previouslyIncluded.add(a3);
              }
            }
          } else {
            // TODO handle higher order bonds and aromaticity.
            //Look for 1_3 atoms that solely attach to a2 (e.g., -C=O-NH2, -C=0-CH3)
            for (Bond b : bondsa2) {
              Atom a3 = b.get1_2(a2);
              if(previouslyIncluded.contains(a3)){
                continue;
              } else if (!a3.isHydrogen() && nonHbonds == 1) {
                int numBondsA3 = a3.getNumBonds();
                if((numBondsA3 == 1) && (a3.isBonded(a2) || a3.isBonded(a))){
                  queueAtoms.add(a3.getIndex());
                  previouslyIncluded.add(a3);
                }
              }
            }
          }
        }
        includeAtoms.addAll(queueAtoms);
        excludeAtoms = "";
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
          int ind = i+1;
          if (!includeAtoms.contains(ind) && excludeAtoms.length() == 0) {
            excludeAtoms = ind;
          } else if (!includeAtoms.contains(ind)) {
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
        if(logger.isLoggable(Level.FINE)){
          logger.fine("A: " + excludeAtomsA);
          logger.fine("B: " + excludeAtomsB);
        }

        // Compare structures in baseFilter and targetFilter.
        ProgressiveAlignmentOfCrystals pac = new ProgressiveAlignmentOfCrystals(baseFilter, targetFilter,
                isSymmetric)
        runningStatistics =
                pac.comparisons(numAU, inflationFactor, matchTol, hitTol, zPrime1, zPrime2, excludeAtomsA, excludeAtomsB,
                        alphaCarbons, includeHydrogen, massWeighted, crystalPriority, strict, saveClusters, save,
                        restart, write, machineLearning, inertia, gyrationComponents, linkage, printSym,
                        lowMemory, createFE, appendSym, pacFilename)
      }
    }else {
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
                      alphaCarbons, includeHydrogen, massWeighted, crystalPriority, strict, saveClusters, save,
                      restart, write, machineLearning, inertia, gyrationComponents, linkage, printSym,
                      lowMemory, createFE, appendSym, pacFilename)
    }
    baseFilter.closeReader();
    targetFilter.closeReader();
    return this
  }
}