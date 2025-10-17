//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.algorithms.commands;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.numerics.Potential;
import ffx.numerics.math.RunningStatistics;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.utils.ProgressiveAlignmentOfCrystals;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;

import static ffx.utilities.StringUtils.parseAtomRanges;
import static java.lang.String.format;
import static org.apache.commons.io.FilenameUtils.concat;
import static org.apache.commons.io.FilenameUtils.getBaseName;
import static org.apache.commons.io.FilenameUtils.getFullPath;

@Command(description = " Determine the RMSD for crystal polymorphs using the Progressive Alignment of Crystals (PAC) algorithm.",
    name = "SuperposeCrystals")
public class SuperposeCrystals extends AlgorithmsCommand {

  @Option(names = {"--ac", "--alchemicalAtoms"}, paramLabel = "", defaultValue = "",
      description = "Atom indices to be excluded from both crystals (e.g. 1-24,32-65). Use if molecular identity and atom ordering are the same for both crystals (otherwise use \"--ac1\" and \"--ac2\".")
  private String excludeAtoms = "";

  @Option(names = {"--ac1", "--alchemicalAtoms1"}, paramLabel = "", defaultValue = "",
      description = "Atom indices to be excluded in the first crystal (e.g. 1-24,32-65).")
  private String excludeAtomsA = "";

  @Option(names = {"--ac2", "--alchemicalAtoms2"}, paramLabel = "", defaultValue = "",
      description = "Atom indices to be excluded in the second crystal (e.g. 1-24,32-65).")
  private String excludeAtomsB = "";

  @Option(names = {"--ca", "--carbonAlphas"}, paramLabel = "false", defaultValue = "false",
      description = "Consider only alpha carbons for proteins.")
  private static boolean alphaCarbons;

  @Option(names = {"--cd", "--createDirectories"}, paramLabel = "false", defaultValue = "false",
      description = "Create subdirectories for free energy simulations.")
  private static boolean createFE;

  @Option(names = {"--as", "--autoSymmetry"}, paramLabel = "false", defaultValue = "false",
      description = "Automatically generate symmetry operators.")
  private static boolean autoSym;

  @Option(names = {"--st", "--symTolerance"}, paramLabel = "0.5", defaultValue = "0.5",
      description = "Add atoms to torsion based symmetries that are beneath this cutoff.")
  private static double symTolerance;

  @Option(names = {"--ws", "--writeSym"}, paramLabel = "false", defaultValue = "false",
      description = "Write the created symmetry operators in the corresponding Properties/Key file.")
  private static boolean writeSym;

  @Option(names = {"--gc", "--gyrationComponents"}, paramLabel = "false", defaultValue = "false",
      description = "Display components for radius of gyration for final clusters.")
  private static boolean gyrationComponents;

  @Option(names = {"--ht", "--hitTolerance"}, paramLabel = "-1.0", defaultValue = "-1.0",
      description = "Sum comparisons that attain a value lower than this tolerance.")
  private double hitTol;

  @Option(names = {"--if", "--inflationFactor"}, paramLabel = "5.0", defaultValue = "5.0",
      description = "Inflation factor used to determine replicates expansion (IF * nAU in replicates).")
  private double inflationFactor;

  @Option(names = {"--ih", "--includeHydrogen"}, paramLabel = "false", defaultValue = "false",
      description = "Include hydrogen atoms.")
  private boolean includeHydrogen;

  @Option(names = {"--in", "--inertia"}, paramLabel = "false", defaultValue = "false",
      description = "Display moments of inertia for final clusters.")
  private static boolean inertia;

  @Option(names = {"-l", "--linkage"}, paramLabel = "1", defaultValue = "1",
      description = "Single (0), Average (1), or Complete (2) coordinate linkage for molecule prioritization.")
  private int linkage;

  @Option(names = {"--lm", "--lowMemory"}, paramLabel = "false", defaultValue = "false",
      description = "Reduce memory usage at the cost of efficiency.")
  private static boolean lowMemory;

  @Option(names = {"--mt", "--moleculeTolerance"}, paramLabel = "0.0015", defaultValue = "0.0015",
      description = "Tolerance to determine if two AUs are different.")
  private double matchTol;

  @Option(names = {"--mw", "--massWeighted"}, paramLabel = "false", defaultValue = "false",
      description = "Use mass-weighted atomic coordinates for alignment.")
  private static boolean massWeighted;

  @Option(names = {"--na", "--numAU"}, paramLabel = "20", defaultValue = "20",
      description = "Number of asymmetric units included to calculate the RMSD.")
  private int numAU;

  @Option(names = {"--pc", "--prioritizeCrystals"}, paramLabel = "0", defaultValue = "0",
      description = "Prioritize crystals based on high density (0), low density (1), or file order (2).")
  private int crystalPriority;

  @Option(names = {"--ps", "--printSymOp"}, paramLabel = "-1.0", defaultValue = "-1.0",
      description = "Print optimal SymOp to align input crystals (print out atom deviations above value).")
  private static double printSym;

  @Option(names = {"-r", "--restart"}, paramLabel = "false", defaultValue = "false",
      description = "Restart from a previously written RMSD matrix (if one exists).")
  private static boolean restart;

  @Option(names = {"-s", "--save"}, paramLabel = "-1.0", defaultValue = "-1.0",
      description = "Save structures less than or equal to this cutoff.")
  private static double save;

  @Option(names = {"--sc", "--saveClusters"}, paramLabel = "0", defaultValue = "0",
      description = "Save files for the superposed crystal clusters (1=PDB, 2=XYZ).")
  private static int saveClusters;

  @Option(names = {"--sm", "--saveMachineLearning"}, paramLabel = "false", defaultValue = "false",
      description = "Final structures for each comparison will be written out with RMSD in a CSV.")
  private static boolean machineLearning;

  @Option(names = {"--th", "--thorough"}, paramLabel = "false", defaultValue = "false",
      description = "More intensive, less efficient version of PAC.")
  private static boolean thorough;

  @Option(names = {"-w", "--write"}, paramLabel = "false", defaultValue = "false",
      description = "Write out the PAC RMSD matrix.")
  private static boolean write;

  @Option(names = {"--zp", "--zPrime"}, paramLabel = "-1", defaultValue = "-1",
      description = "Z'' for both crystals (assumes same value).")
  private int zPrime;

  @Option(names = {"--zp1", "--zPrime1"}, paramLabel = "-1", defaultValue = "-1",
      description = "Z'' for crystal 1 (default will try to autodetect).")
  private int zPrime1;

  @Option(names = {"--zp2", "--zPrime2"}, paramLabel = "-1", defaultValue = "-1",
      description = "Z'' for crystal 2 (default will try to autodetect).")
  private int zPrime2;

  @Parameters(arity = "1..2", paramLabel = "files",
      description = "Atomic coordinate file(s) to compare in XYZ format.")
  List<String> filenames = null;

  private final List<Atom> previouslyIncludedBase = new ArrayList<>();
  private final List<Integer> queueIndices = new ArrayList<>();
  public RunningStatistics runningStatistics;
  public final StringBuilder symOpsA = new StringBuilder();
  public final StringBuilder symOpsB = new StringBuilder();

  public SuperposeCrystals() {
    super();
  }

  public SuperposeCrystals(FFXBinding binding) {
    super(binding);
  }

  public SuperposeCrystals(String[] args) {
    super(args);
  }

  @Override
  public SuperposeCrystals run() {
    System.setProperty("vdwterm", "false");

    if (!init()) {
      return this;
    }

    if (filenames == null) {
      logger.info(helpString());
      return this;
    }

    algorithmFunctions.openAll(filenames.get(0));
    SystemFilter baseFilter = algorithmFunctions.getFilter();
    SystemFilter targetFilter;

    boolean isSymmetric = false;
    int numFiles = filenames.size();
    if (numFiles == 1) {
      logger.info(
          "\n PAC will be applied between all pairs of conformations within the supplied file.\n");
      isSymmetric = true;
      algorithmFunctions.openAll(filenames.get(0));
      targetFilter = algorithmFunctions.getFilter();
    } else {
      logger.info(
          "\n PAC will compare all conformations in the first file to all those in the second file.\n");
      algorithmFunctions.openAll(filenames.get(1));
      targetFilter = algorithmFunctions.getFilter();
      String filenameA = getBaseName(filenames.get(0));
      String filenameB = getBaseName(filenames.get(1));
      symOpsA.append(format("\n Add the following symop to properties/key file of %s (Apply on %s to generate %s):\nsymop", filenameA, filenameB, filenameA));
      symOpsB.append(format("\n\n Add the following symop to properties/key file of %s (Apply on %s to generate %s):\nsymop", filenameB, filenameA, filenameB));
    }

    String filename = filenames.get(0);
    String pacFilename = concat(getFullPath(filename), getBaseName(filename) + ".dst");

    if (zPrime > 0 && zPrime % 1 == 0) {
      zPrime1 = zPrime;
      zPrime2 = zPrime;
    }

    if (printSym < 0 && writeSym || autoSym) {
      logger.info("\n Printing distance for atoms greater than large distance of 10.0 Å");
      printSym = 10.0;
    }

    if (createFE) {
      lowMemory = true;
    }

    if (autoSym) {
      boolean loggerLevelChanged = false;
      if (logger.isLoggable(Level.INFO) && !logger.isLoggable(Level.FINE)) {
        loggerLevelChanged = true;
        logger.setLevel(Level.WARNING);
      }

      final Atom[] atomsBase = baseFilter.getActiveMolecularSystem().getAtomArray();
      final Atom[] atomsTarget = targetFilter.getActiveMolecularSystem().getAtomArray();
      final int nAtomsBase = atomsBase.length;
      final int nAtomsTarget = atomsTarget.length;

      List<List<Atom>> smallAtomGroups = new ArrayList<>();
      List<Double> smallGroupsRMSD = new ArrayList<>();

      List<Atom> mapBaseAtoms = new ArrayList<>();
      List<Integer> mapBaseIndices = new ArrayList<>();
      List<Atom> mapTargetAtoms = new ArrayList<>();
      List<Integer> mapTargetIndices = new ArrayList<>();

      if (!excludeAtoms.isEmpty()) {
        List<Integer> alchAtoms = parseAtomRanges("AutoSym", excludeAtoms, nAtomsBase);
        for (int i = 0; i < nAtomsBase - 1; i++) {
          if (!alchAtoms.contains(i + 1)) {
            logger.warning(format(" Include Atoms for both:\n %s\n %s.", atomsBase[i].toString(), atomsTarget[i].toString()));
            mapBaseIndices.add(i + 1);
            mapBaseAtoms.add(atomsBase[i]);
            mapTargetIndices.add(i + 1);
            mapTargetAtoms.add(atomsTarget[i]);
          } else {
            previouslyIncludedBase.add(atomsBase[i]);
          }
        }
      } else if (!excludeAtomsA.isEmpty() || !excludeAtomsB.isEmpty()) {
        logger.warning(format(" Independent alchemical atoms has not been implemented."));
      }

      assert (mapBaseIndices.size() == mapBaseAtoms.size());
      assert (mapTargetIndices.size() == mapTargetAtoms.size());
      if (mapBaseAtoms.size() != mapTargetAtoms.size()) {
        logger.warning(format(" Error in atom selection. Base size of %3d does not match Target size of %3d.", mapBaseAtoms.size(), mapTargetAtoms.size()));
        if (logger.isLoggable(Level.FINE)) {
          StringBuilder sb = new StringBuilder();
          for (Integer ind : mapBaseIndices) {
            sb.append(ind).append(",");
          }
          sb.append("\n");
          logger.fine(" Base indices:" + sb.toString());
          sb.setLength(0);
          for (Integer ind : mapTargetIndices) {
            sb.append(ind).append(",");
          }
          sb.append("\n");
          logger.fine(" Target indices:" + sb.toString());
        }
        return this;
      }

      for (Atom a : atomsBase) {
        if (previouslyIncludedBase.contains(a)) {
          continue;
        }
        queueIndices.clear();
        addAtomIndex(a);

        ArrayList<Atom> queueAtoms = new ArrayList<>();
        for (Integer index : queueIndices) {
          queueAtoms.add(atomsBase[index - 1]);
        }

        if (queueAtoms.size() < 3) {
          boolean groupFound = false;
          for (Atom atom : queueAtoms) {
            Bond[] bonds = atom.getBonds().toArray(new Bond[0]);
            for (Bond bond : bonds) {
              Atom atomGrouped = bond.get1_2(atom);
              for (List<Atom> atomList : smallAtomGroups) {
                if (atomList.contains(atomGrouped)) {
                  groupFound = true;
                  for (Atom atomNeedingGroup : queueAtoms) {
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
          smallAtomGroups.add(queueAtoms);
        }
      }

      for (int i = 0; i < nAtomsBase; i++) {
        Atom atom = atomsBase[i];
        boolean found = false;
        for (List<Atom> list : smallAtomGroups) {
          if (list.contains(atom)) {
            found = true;
            break;
          }
        }
        if (!found) {
          boolean groupFound = false;
          Bond[] bonds = atom.getBonds().toArray(new Bond[0]);
          for (Bond bond : bonds) {
            Atom atomGrouped = bond.get1_2(atom);
            for (List<Atom> atomList : smallAtomGroups) {
              if (atomList.contains(atomGrouped)) {
                groupFound = true;
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

      StringBuilder sb = new StringBuilder(" Small Groups:\n");
      for (List<Atom> group : smallAtomGroups) {
        ArrayList<Integer> includeIndices = new ArrayList<>();
        sb.append(" New Group: ");
        for (Atom atom : group) {
          sb.append(atom.getIndex() + ", ");
          includeIndices.add(atom.getIndex());
        }
        sb.append("\n");

        excludeAtoms = "";

        for (int i = 0; i < nAtomsBase; i++) {
          int ind = i + 1;
          if (!includeIndices.contains(ind) && excludeAtoms.length() == 0) {
            excludeAtoms = String.valueOf(ind);
          } else if (!includeIndices.contains(ind)) {
            excludeAtoms += "," + ind;
          }
        }

        if (excludeAtoms == null || excludeAtoms.isEmpty()) {
          continue;
        } else {
          excludeAtomsA = excludeAtoms;
          excludeAtomsB = excludeAtoms;
        }
        if (logger.isLoggable(Level.FINE)) {
          logger.fine("A: " + excludeAtomsA);
          logger.fine("B: " + excludeAtomsB);
        }

        ProgressiveAlignmentOfCrystals pac = new ProgressiveAlignmentOfCrystals(baseFilter, targetFilter,
            isSymmetric);
        runningStatistics =
            pac.comparisons(numAU, inflationFactor, matchTol, hitTol, zPrime1, zPrime2, excludeAtomsA, excludeAtomsB,
                alphaCarbons, includeHydrogen, massWeighted, crystalPriority, thorough, saveClusters, save,
                restart, write, machineLearning, inertia, gyrationComponents, linkage, printSym,
                lowMemory, createFE, false, pacFilename, new StringBuilder(""), new StringBuilder(""));
        double min = runningStatistics.getMin();
        smallGroupsRMSD.add(min);
      }

      logger.info(sb.toString());
      int numGroups = smallAtomGroups.size();
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Number of minimal groups: %3d", numGroups));
      }

      int[][] smallIndexGroups = new int[numGroups][];
      for (int i = 0; i < numGroups; i++) {
        List<Atom> currentAtoms = smallAtomGroups.get(i);
        int currentSize = currentAtoms.size();
        smallIndexGroups[i] = new int[currentSize];
        for (int j = 0; j < currentSize; j++) {
          smallIndexGroups[i][j] = currentAtoms.get(j).getIndex();
        }
      }

      ArrayList<ArrayList<Integer>> finalGroups = new ArrayList<>();
      ArrayList<Integer> acceptedAtoms = new ArrayList<>();
      ArrayList<Integer> queuedIndices = new ArrayList<>();
      ArrayList<Integer> usedGroups = new ArrayList<>();

      for (int i = 0; i < numGroups; i++) {
        List<Atom> currentAtomGroup = smallAtomGroups.get(i);
        int[] currentIndexGroup = smallIndexGroups[i];
        int numAtomsInGroup = currentAtomGroup.size();
        acceptedAtoms.clear();
        queuedIndices.clear();
        if (logger.isLoggable(Level.FINE)) {
          logger.fine(format(" Minimal group number: %3d Size: %3d RMSD: %9.3f", i, numAtomsInGroup, smallGroupsRMSD.get(i)));
        }
        if (!usedGroups.contains(i)) {
          for (int index : currentIndexGroup) {
            acceptedAtoms.add(index);
          }
          usedGroups.add(i);

          if (smallGroupsRMSD.get(i) > symTolerance) {
            if (logger.isLoggable(Level.FINE)) {
              logger.fine(format(" Run (group) %3d added outright to final.", i));
            }
          } else {
            for (int j = 0; j < acceptedAtoms.size(); j++) {
              Atom a1 = atomsBase[acceptedAtoms.get(j) - 1];
              if (logger.isLoggable(Level.FINE)) {
                logger.fine(format(" Current group: %3d Atom in Group: %3d of %3d (was %3d==%3d) Atom Index: %3d", i, j, acceptedAtoms.size(), numAtomsInGroup, currentIndexGroup.length, a1.getIndex()));
              }
              Bond[] b1 = a1.getBonds().toArray(new Bond[0]);
              for (Bond b : b1) {
                Atom a2 = b.get1_2(a1);
                int a2Index = a2.getIndex();
                boolean currentIndexGroupContains = false;
                for (int idx : currentIndexGroup) {
                  if (idx == a2Index) {
                    currentIndexGroupContains = true;
                    break;
                  }
                }
                if (!currentIndexGroupContains) {
                  for (int k = 0; k < numGroups; k++) {
                    int[] currentIndexGroup2 = smallIndexGroups[k];
                    boolean currentIndexGroup2Contains = false;
                    for (int idx : currentIndexGroup2) {
                      if (idx == a2Index) {
                        currentIndexGroup2Contains = true;
                        break;
                      }
                    }
                    if (!usedGroups.contains(k) && smallGroupsRMSD.get(k) < symTolerance && currentIndexGroup2Contains) {
                      for (int ind : currentIndexGroup2) {
                        queuedIndices.add(ind);
                      }

                      excludeAtoms = "";
                      for (int l = 0; l < nAtomsBase; l++) {
                        int ind = l + 1;
                        boolean exclude = !acceptedAtoms.contains(ind) && !queuedIndices.contains(ind);
                        if (exclude) {
                          if (excludeAtoms.length() == 0) {
                            excludeAtoms = String.valueOf(ind);
                          } else {
                            excludeAtoms += "," + ind;
                          }
                        }
                      }

                      excludeAtomsA = excludeAtoms;
                      excludeAtomsB = excludeAtoms;

                      if (logger.isLoggable(Level.FINE)) {
                        logger.fine(" Excluded atoms A: " + excludeAtomsA);
                        logger.fine(" Excluded atoms B: " + excludeAtomsB);
                      }

                      ProgressiveAlignmentOfCrystals pac = new ProgressiveAlignmentOfCrystals(baseFilter, targetFilter,
                          isSymmetric);
                      runningStatistics =
                          pac.comparisons(numAU, inflationFactor, matchTol, hitTol, zPrime1, zPrime2, excludeAtomsA, excludeAtomsB,
                              alphaCarbons, includeHydrogen, massWeighted, crystalPriority, thorough, saveClusters, save,
                              restart, write, machineLearning, inertia, gyrationComponents, linkage, printSym,
                              lowMemory, createFE, false, pacFilename, new StringBuilder(""), new StringBuilder(""));
                      double min = runningStatistics.getMin();
                      if (min > symTolerance) {
                        if (logger.isLoggable(Level.FINE)) {
                          logger.fine(format(" Run %3d Group %3d failed with RMSD %9.3f > tolerance (%9.3f).", i, k, min, symTolerance));
                        }
                        queuedIndices.clear();
                      } else {
                        for (Integer value : queuedIndices) {
                          if (!acceptedAtoms.contains(value)) {
                            acceptedAtoms.add(value);
                          }
                        }
                        if (logger.isLoggable(Level.FINE)) {
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
          ArrayList<Integer> clone = new ArrayList<>(acceptedAtoms.size());
          clone.addAll(acceptedAtoms);
          finalGroups.add(clone);
          if (logger.isLoggable(Level.FINE)) {
            logger.fine(format(" Run %3d accepted atoms (size: %3d) added to finalGroups (size: %3d)", i, acceptedAtoms.size(), finalGroups.size()));
          }
        }
      }

      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Final Num Groups: %3d", finalGroups.size()));
      }

      if (loggerLevelChanged) {
        logger.setLevel(Level.INFO);
      }

      for (ArrayList<Integer> group : finalGroups) {
        int groupSize = group.size();
        if (logger.isLoggable(Level.FINE)) {
          logger.fine(format(" Final Size: %3d", groupSize));
        }
        if (groupSize <= 0) {
          logger.warning(" Final group of size zero was encountered. Check selection logic.");
          continue;
        }
        queuedIndices.clear();
        for (Integer ind : group) {
          queuedIndices.add(ind);
        }
        excludeAtoms = "";
        for (int l = 0; l < nAtomsBase; l++) {
          int ind = l + 1;
          if (!queuedIndices.contains(ind) && excludeAtoms.length() == 0) {
            excludeAtoms = String.valueOf(ind);
          } else if (!queuedIndices.contains(ind)) {
            excludeAtoms += "," + ind;
          }
        }
        excludeAtomsA = excludeAtoms;
        excludeAtomsB = excludeAtoms;
        if (logger.isLoggable(Level.FINE)) {
          logger.fine("A: " + excludeAtomsA);
          logger.fine("B: " + excludeAtomsB);
        }

        ProgressiveAlignmentOfCrystals pac = new ProgressiveAlignmentOfCrystals(baseFilter, targetFilter,
            isSymmetric);
        runningStatistics =
            pac.comparisons(numAU, inflationFactor, matchTol, hitTol, zPrime1, zPrime2, excludeAtomsA, excludeAtomsB,
                alphaCarbons, includeHydrogen, massWeighted, crystalPriority, thorough, saveClusters, save,
                restart, write, machineLearning, inertia, gyrationComponents, linkage, printSym,
                lowMemory, createFE, writeSym, pacFilename, symOpsA, symOpsB);
      }
    } else {
      if (excludeAtoms != null && !excludeAtoms.isEmpty()) {
        excludeAtomsA = excludeAtoms;
        excludeAtomsB = excludeAtoms;
      }
      ProgressiveAlignmentOfCrystals pac = new ProgressiveAlignmentOfCrystals(baseFilter, targetFilter,
          isSymmetric);
      runningStatistics =
          pac.comparisons(numAU, inflationFactor, matchTol, hitTol, zPrime1, zPrime2, excludeAtomsA, excludeAtomsB,
              alphaCarbons, includeHydrogen, massWeighted, crystalPriority, thorough, saveClusters, save,
              restart, write, machineLearning, inertia, gyrationComponents, linkage, printSym,
              lowMemory, createFE, writeSym, pacFilename, symOpsA, symOpsB);
    }

    if (printSym >= 0.0) {
      int finalCharA = symOpsA.length() - 2;
      if (symOpsA.substring(finalCharA).equals("\\\n")) {
        symOpsA.setLength(finalCharA);
      }
      int finalCharB = symOpsB.length() - 2;
      if (symOpsB.substring(finalCharB).equals("\\\n")) {
        symOpsB.setLength(finalCharB);
      }
      if (symOpsA.indexOf("\\") >= 0 && symOpsB.indexOf("\\") >= 0) {
        logger.info(symOpsA.toString() + symOpsB.toString());
      } else if (logger.isLoggable(Level.FINE)) {
        logger.fine(" Following generated SymOp strings did not contain symmetry operators." + symOpsA.toString() + symOpsB.toString());
      }
    }

    baseFilter.closeReader();
    targetFilter.closeReader();
    return this;
  }

  private void addAtomIndex(Atom aAdded) {
    if (!previouslyIncludedBase.contains(aAdded)) {
      queueIndices.add(aAdded.getIndex());
      previouslyIncludedBase.add(aAdded);
      int numBondsAdded = aAdded.getNumBonds();
      int numHAdded = aAdded.getNumberOfBondedHydrogen();
      int nonHbonds = numBondsAdded - numHAdded;

      Bond[] bonds4Added = aAdded.getBonds().toArray(new Bond[0]);
      for (Bond b : bonds4Added) {
        Atom aBound2Added = b.get1_2(aAdded);
        if (nonHbonds == 1) {
          addAtomIndex(aBound2Added);
          for (Bond bAdded : bonds4Added) {
            Atom a3 = bAdded.get1_2(aAdded);
            if (a3.isHydrogen()) {
              addAtomIndex(a3);
            }
          }
        }

        int numBondsB2A = aBound2Added.getNumBonds();
        int numHB2A = aBound2Added.getNumberOfBondedHydrogen();
        int nonHB2A = numBondsB2A - numHB2A;
        Bond[] bondsB2A = aBound2Added.getBonds().toArray(new Bond[0]);
        if (nonHB2A == 1) {
          addAtomIndex(aBound2Added);
          for (Bond b2a : bondsB2A) {
            Atom aBond2B2A = b2a.get1_2(aBound2Added);
            if (aBond2B2A.isHydrogen()) {
              addAtomIndex(aBond2B2A);
            }
          }
        } else if (numBondsAdded == 1 || numBondsB2A == 1) {
          addAtomIndex(aBound2Added);
        } else if (nonHB2A == 2) {
          addAtomIndex(aBound2Added);
          for (Bond b2a : bondsB2A) {
            Atom aBond2B2A = b2a.get1_2(aBound2Added);
            if (aBond2B2A.isHydrogen()) {
              addAtomIndex(aBond2B2A);
            }
          }
        } else {
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

  @Override
  public List<Potential> getPotentials() {
    return Collections.emptyList();
  }
}
