// ******************************************************************************
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
// ******************************************************************************
package ffx.algorithms.cli;

import static java.lang.String.format;

import ffx.algorithms.optimize.RotamerOptimization;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.*;
import ffx.potential.cli.PotentialScript;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.TitrationUtils;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.logging.Logger;

import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.configuration2.Configuration;
import org.apache.commons.io.FilenameUtils;
import picocli.CommandLine.ArgGroup;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that use a many-body expansion for global
 * optimization.
 *
 * @author Michael J. Schnieders
 * @author Mallory R. Tollefson
 * @since 1.0
 */
public class ManyBodyOptions{

  private static final Logger logger = Logger.getLogger(ManyBodyOptions.class.getName());

  /**
   * The ArgGroup keeps the ManyBodyOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Many-Body Optimization Options%n", validate = false)
  public ManyBodyOptionGroup group = new ManyBodyOptionGroup();

  /**
   * The ArgGroup keeps the ManyBodyBoxOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Many-Body Box Optimization Options%n", validate = false)
  public ManyBodyBoxOptionGroup boxGroup = new ManyBodyBoxOptionGroup();

  /**
   * The ArgGroup keeps the ManyBodyWindowOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Many-Body Window Optimization Options%n", validate = false)
  public ManyBodyWindowOptionGroup windowGroup = new ManyBodyWindowOptionGroup();

  /**
   * The ArgGroup keeps the ManyBodyWindowOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Many-Body Energy Expansion and Cut-off Options%n", validate = false)
  public ManyBodyEnergyOptionGroup energyGroup = new ManyBodyEnergyOptionGroup();

  /**
   * The ArgGroup keeps the ManyBodyResidueOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Many-Body Residue Selection Options%n", validate = false)
  public ManyBodyResidueOptionGroup residueGroup = new ManyBodyResidueOptionGroup();

  private RotamerOptimization rotamerOptimization;
  private RotamerLibrary rotamerLibrary;
  private int allStartResID;
  private int boxStart;
  private int boxEnd;
  private int[] numXYZBoxes;
  private int forceResiduesStart;
  private int forceResiduesEnd;

  /**
   * Gets the algorithm number.
   *
   * @return Integer representing the algorithm being run (i.e. global, box optimization, window,
   *     etc.)
   */
  public int getAlgorithmNumber() {
    return group.algorithm;
  }

  public boolean getUsingOriginalCoordinates() {
    return !group.noOriginal;
  }

  /**
   * initRotamerOptimization.
   *
   * @param rotamerOptimization a {@link RotamerOptimization} object.
   * @param activeAssembly a {@link ffx.potential.MolecularAssembly} object.
   */
  public void initRotamerOptimization(
          RotamerOptimization rotamerOptimization, MolecularAssembly activeAssembly, String fileName) {
    this.rotamerOptimization = rotamerOptimization;
    logger.info("Init rotopt");
    boolean useOrigCoordsRotamer = !group.noOriginal;
    if (group.decompose) {
      useOrigCoordsRotamer = true;
    }

    String rotamerFileName = activeAssembly.getFile().getName();
    rotamerFileName = FilenameUtils.removeExtension(rotamerFileName);
    rotamerFileName = rotamerFileName + ".rot";
    File rotFile = new File(rotamerFileName);

    if (rotFile.exists()) {
      logger.info(" EXPERIMENTAL: Using rotamer file " + rotamerFileName);
      rotamerLibrary = new RotamerLibrary(RotamerLibrary.ProteinLibrary.None, false);
      try {
        RotamerLibrary.readRotFile(rotFile, activeAssembly);
      } catch (IOException iox) {
        logger.severe(
            format(
                " Exception in parsing rotamer file: %s\n%s",
                iox, Utilities.stackTraceToString(iox)));
      }
    } else {
      rotamerLibrary =
          new RotamerLibrary(
              RotamerLibrary.ProteinLibrary.intToProteinLibrary(group.library),
              useOrigCoordsRotamer);
    }

    rotamerOptimization.setRotamerLibrary(rotamerLibrary);
    rotamerOptimization.setSingletonClashThreshold(energyGroup.clashThreshold);
    rotamerOptimization.setPairClashThreshold(energyGroup.pairClashThreshold);
    rotamerOptimization.setDecomposeOriginal(group.decompose);

    if (group.algorithm == 0) {
      setAlgorithm(activeAssembly);
    }
    setSelection();
    setForcedResidue();
    setResidues(activeAssembly, fileName);
    setRotOptProperties();
  }

  public void setOriginalCoordinates(boolean useOrig) {
    group.noOriginal = !useOrig;
  }

  /**
   * setResidues.
   *
   * @param activeAssembly a {@link ffx.potential.MolecularAssembly} object.
   */
  public void setResidues(MolecularAssembly activeAssembly, String fileName) {
    List<String> resList = new ArrayList<>();
    addListResidues(resList);

    TitrationUtils titrationUtils;
    if (group.titrationPH > 0.0 && group.titrationPH <= 14.0) {
      logger.info(format(" Turning on ASP, GLU, LYS and HIS titration rotamers at pH %5.2f", group.titrationPH));

      titrationUtils = new TitrationUtils(activeAssembly.getForceField());
      titrationUtils.setRotamerPhBias(298.15, group.titrationPH);
      List<Residue> residues = activeAssembly.getResidueList();

      for (Residue residue : residues) {
        residue.setTitrationUtils(titrationUtils);
      }

    }

    int counter = 1;
    if (group.algorithm != 5) {
      if (allStartResID > 0) {
        List<Residue> residueList = getResidues(activeAssembly);
        Polymer[] polymers = activeAssembly.getChains();
        for (Polymer polymer : polymers) {
          List<Residue> residues = polymer.getResidues();
          for (Residue residue : residues) {
            Rotamer[] rotamers = residue.setRotamers(rotamerLibrary);
            if (rotamers != null) {
              int nrot = rotamers.length;
              if (nrot == 1) {
                RotamerLibrary.applyRotamer(residue, rotamers[0]);
              }
              if (counter >= allStartResID) {
                residueList.add(residue);
              }
            } else if (!group.forceResidues.equalsIgnoreCase("none")) {
              if (counter >= allStartResID
                  && counter >= forceResiduesStart
                  && counter <= forceResiduesEnd) {
                residueList.add(residue);
              }
            }
            counter++;
          }
        }

        rotamerOptimization.setResidues(residueList);

      } else if (!residueGroup.listResidues.equalsIgnoreCase("none")) {
        List<Residue> residueList = getResidues(activeAssembly);
        Polymer[] polymers = activeAssembly.getChains();
        int n = 0;
        for (String s : resList) {
          Character chainID = s.charAt(0);
          int i = Integer.parseInt(s.substring(1));
          for (Polymer p : polymers) {
            if (p.getChainID() == chainID) {
              List<Residue> rs = p.getResidues();
              for (Residue r : rs) {
                if (r.getResidueNumber() == i) {
                  residueList.add(r);
                  Rotamer[] rotamers = r.setRotamers(rotamerLibrary);
                  if (rotamers != null) {
                    n++;
                  }
                }
              }
            }
          }
        }
        rotamerOptimization.setResiduesIgnoreNull(residueList);
        if (residueList == null) {
          return;
        }
      } else if (!residueGroup.chain.equalsIgnoreCase("-1")) {
        rotamerOptimization.setResidues(residueGroup.chain, residueGroup.start, residueGroup.finish);
      } else {
        rotamerOptimization.setResidues(residueGroup.start, residueGroup.finish);
      }
    } else {
      List<Residue> residueList = new ArrayList<>();
      Polymer[] polymers = activeAssembly.getChains();

      CompositeConfiguration properties = activeAssembly.getProperties();
      boolean ignoreNA = properties.getBoolean("ignoreNA", false);

      if (!residueGroup.listResidues.equalsIgnoreCase("none")) {
        int n = 0;
        for (String s : resList) {
          Character chainID = s.charAt(0);
          int i = Integer.parseInt(s.substring(1));
          for (Polymer p : polymers) {
            if (p.getChainID() == chainID) {
              List<Residue> rs = p.getResidues();
              for (Residue r : rs) {
                if (ignoreNA && r.getResidueType() == Residue.ResidueType.NA) {
                  continue;
                }
                if (r.getResidueNumber() == i) {
                  residueList.add(r);
                  Rotamer[] rotamers = r.setRotamers(rotamerLibrary);
                  if (rotamers != null) {
                    n++;
                  }
                }
              }
            }
          }
        }
        rotamerOptimization.setResiduesIgnoreNull(residueList);
        if (n < 1) {
          return;
        }
      } else {
        for (Polymer p : polymers) {
          List<Residue> rs = p.getResidues();
          for (Residue r : rs) {
            if (ignoreNA && r.getResidueType() == Residue.ResidueType.NA) {
              continue;
            }
            Rotamer[] rotamers = r.setRotamers(rotamerLibrary);
            if (rotamers != null) {
              int nrot = rotamers.length;
              if (nrot == 1) {
                RotamerLibrary.applyRotamer(r, rotamers[0]);
              } else if (nrot > 1) {
                residueList.add(r);
              }
            }
            counter++;
          }
        }

        //            boolean ignoreNA = false;
        //            String ignoreNAProp = System.getProperty("ignoreNA");
        //            if (ignoreNAProp != null && ignoreNAProp.equalsIgnoreCase("true")) {
        //                ignoreNA = true;
        //            }
        //            List<Residue> residueList = new ArrayList<>();
        //            Polymer[] polymers = activeAssembly.getChains();
        //            int nPolymers = polymers.length;
        //            for (int p = 0; p < nPolymers; p++) {
        //                Polymer polymer = polymers[p];
        //                List<Residue> residues = polymer.getResidues();
        //
        //                System.out.print("\nresidues:\n");
        //
        //                int nResidues = residues.size();
        //                for (int i = 0; i < nResidues; i++) {
        //                    Residue residue = residues.get(i);
        //
        //                    System.out.print(residue+"\n");
        //
        //                    if (ignoreNA && residue.getResidueType() == ResidueType.NA) {
        //                        continue;
        //                    }
        //                    Rotamer[] rotamers = residue.getRotamers(rotamerLibrary);
        //                    if (rotamers != null) {
        //                        int nrot = rotamers.length;
        //                        if (nrot == 1) {
        //                            RotamerLibrary.applyRotamer(residue, rotamers[0]);
        //                        } else if (nrot > 1) {
        //                            residueList.add(residue);
        //                        }
        //                    }
        //                    counter++;
        //                }
      }
      rotamerOptimization.setResidues(residueList);
      rotamerOptimization.setBoxStart(boxStart);
      if (residueGroup.finish > 0) {
        rotamerOptimization.setBoxEnd(boxEnd);
      }
    }

  }

  public List<Residue> getResidues(MolecularAssembly activeAssembly){
    List<Residue> residueList = new ArrayList<>();
    int counter = 0;
    List<Residue> residues = activeAssembly.getResidueList();
    if (residueGroup.all > -1) {
      counter = residueGroup.all;
      for(Residue residue : residues){
        if (residue.getResidueNumber() == counter){
          residueList.add(residue);
          counter += 1;
        }
      }
    } else if (residueGroup.start > -1) {
      counter = residueGroup.start;
      for(Residue residue : residues){
        if (counter == residueGroup.finish + 1){
          break;
        }else if (residue.getResidueNumber() == counter){
          residueList.add(residue);
          counter += 1;

        }
      }
    } else {
      residueList = residues;
    }
    return residueList;
  }

  /**
   * This method sets the algorithm by default. If no parameters are given, the default algorithm
   * value 0. When the default algorithm is 0, a specific algorithm number (1-5) needs to be
   * assigned. This method assigns the default algorithm based on variation in input parameters
   * (start, finish, all, chain, etc.) and then determines how many amino acids are to be optimized.
   * If more than 100 amino acids are to be optimized, the algorithm is set to use box optimization.
   * If fewer than 100 amino acids are to be optimized, the algorithm is set to use global
   * optimization.
   *
   * @param activeAssembly The protein to be optimized.
   * @return The value that the algorithm should be set to since no 1-5 value was assigned by the
   *     user. Either set to 2-global or 5-box optimization.
   */
  private int setAlgorithm(MolecularAssembly activeAssembly) {

    setStartAndEndDefault();

    if (allStartResID > 0) {
      Polymer[] polymers = activeAssembly.getChains();
      int nResidues = 0;
      for (Polymer polymer : polymers) {
        List<Residue> residues = polymer.getResidues();
        nResidues = residues.size() + nResidues;
      }
      if (nResidues > 100) {
        group.algorithm = 5;
      } else {
        group.algorithm = 2;
      }
    } else if (!residueGroup.listResidues.equalsIgnoreCase("none")) {
      List<String> resList = new ArrayList<>();
      addListResidues(resList);
      if (resList.size() > 100) {
        group.algorithm = 5;
      } else {
        group.algorithm = 2;
      }
    } else if (!residueGroup.chain.equalsIgnoreCase("-1")) {
      group.algorithm = startFinishDifference();
    } else if (residueGroup.start > 0 && residueGroup.finish > 0) {
      group.algorithm = startFinishDifference();
    }
    return group.algorithm;
  }

  public double getApproximate() {
    return rotamerOptimization.getApproximate();
  }

  /**
   * Gets the restart file created during rotamer optimization.
   *
   * @return The restart file.
   */
  public File getRestartFile() {
    return rotamerOptimization.getRestartFile();
  }

  public RotamerLibrary getRotamerLibrary() {
    return rotamerLibrary;
  }

  /**
   * This method calculates the difference between the start and finish variables and it returns the
   * algorithm that should be used by default. If more than 100 side-chains are to be optimized, the
   * default algorithm is box optimization (algorithm = 5). If fewer than 100 residues are to be
   * optimized, the default algorithm is global optimization (algorithm = 2).
   *
   * @return The default algorithm number.
   */
  private int startFinishDifference() {
    if (residueGroup.finish - residueGroup.start > 100) {
      return 5;
    } else {
      return 2;
    }
  }

  /**
   * This method sets the start and finish points by default. If no parameters are given as input,
   * the default behavior is to begin the optimization at the first available amino acid and end the
   * optimization at the last amino acid. The default algorithm (i.e. box versus global) depends on
   * how many amino acids will be optimized (general cutoff of 100 amino acids).
   */
  private void setStartAndEndDefault() {
    if (residueGroup.start < 0 && residueGroup.finish < 0 && residueGroup.all < 0) {
      if (!residueGroup.listResidues.equalsIgnoreCase("none")) {
        allStartResID = -1;
        boxStart = 0;
        boxEnd = -1;
      } else {
        allStartResID = 1;
        boxStart = residueGroup.start - 1;
        boxEnd = residueGroup.finish - 1;
      }
    } else {
      allStartResID = residueGroup.all;
      boxStart = residueGroup.start - 1;
      boxEnd = residueGroup.finish - 1;
    }
  }

  /** Set allStartResID, boxStart and boxEnd */
  private void setSelection() {
    // Chain, Residue and/or Box selections.
    // Internal machinery indexed 0 to (n-1)
    setStartAndEndDefault();

    if (group.algorithm != 5) {
      // Not Box optimization.
      if (allStartResID < 1 && residueGroup.listResidues.equalsIgnoreCase("none")) {
        if (residueGroup.finish < residueGroup.start || residueGroup.start < 0
            || residueGroup.finish < 0) {
          logger.warning(" FFX shutting down: no residues specified for optimization.");
          return;
        }
      }
    } else {
      // Box optimization.
      if (allStartResID > 0) {
        // Internal machinery indexed 0 to (n-1)
        boxStart = allStartResID - 1;
      } else {
        if (boxStart < 0 || (boxEnd > -1 && boxEnd < boxStart)) {
          logger.warning(
              " FFX shutting down: Invalid input for box selection: index begins at 1 (or start must be less than finish).");
          return;
        }
      }
    }

    // Box optimization options.
    numXYZBoxes = new int[3];
    if (group.algorithm == 5) {
      String input = boxGroup.numBoxes;
      Scanner boxNumInput = new java.util.Scanner(input);
      boxNumInput.useDelimiter(",");
      int inputLoopCounter = 0;
      numXYZBoxes[0] = 3; // Default
      while (inputLoopCounter < numXYZBoxes.length) {
        if (boxNumInput.hasNextInt()) {
          numXYZBoxes[inputLoopCounter] = boxNumInput.nextInt();
          inputLoopCounter++;
        } else if (boxNumInput.hasNextDouble()) {
          numXYZBoxes[inputLoopCounter] = (int) Math.floor(boxNumInput.nextDouble());
          inputLoopCounter++;
          logger.info(" Double input to nB truncated to integer.");
        } else if (boxNumInput.hasNext()) {
          logger.info(" Non-numeric input to nB discarded");
          boxNumInput.next();
        } else {
          logger.info(
              " Insufficient input to nB. Non-input values assumed either equal to X or default to 3");
          break;
        }
      }
      boxNumInput.close();
      for (int i = inputLoopCounter; i < numXYZBoxes.length; i++) {
        numXYZBoxes[i] = numXYZBoxes[0];
      }
      for (int i = 0; i < numXYZBoxes.length; i++) {
        if (numXYZBoxes[i] == 0) {
          numXYZBoxes[i] = 3;
          logger.info(" Input of zero to nB reset to default of three.");
        } else if (numXYZBoxes[i] < 0) {
          numXYZBoxes[i] = -1 * numXYZBoxes[i];
          logger.info(" Input of negative number to nB reset to positive number");
        }
      }
    }
  }

  /** setForcedResidue. */
  private void setForcedResidue() {
    // Force residues.
    forceResiduesStart = -1;
    forceResiduesEnd = -1;

    List<String> resList = new ArrayList<>();
    if (!residueGroup.listResidues.equalsIgnoreCase("none")) {
      String[] tok = residueGroup.listResidues.split(",");
      for (String t : tok) {
        logger.info(" Adding " + t);
        resList.add(t);
      }
    }

    // Evaluate forced residues for the sliding window algorithm
    if (group.algorithm == 4 && !group.forceResidues.equalsIgnoreCase("-1,-1")) {
      String input = group.forceResidues;
      Scanner frScan = new Scanner(input);
      frScan.useDelimiter(",");
      try {
        if (!frScan.hasNextInt()) {
          frScan.next(); // Discards extra input to indicate a negative value of frStart.
        }
        forceResiduesStart = frScan.nextInt();
        forceResiduesEnd = frScan.nextInt();
      } catch (Exception ex) {
        logger.severe(
            format(
                " FFX shutting down: input to -fR could not be parsed as a pair of integers: %s",
                input));
      }
      if (forceResiduesStart > forceResiduesEnd) {
        logger.info(" Start of range higher than ending: start flipped with end.");
        int temp = forceResiduesStart;
        forceResiduesStart = forceResiduesEnd;
        forceResiduesEnd = temp;
      }
      if (forceResiduesEnd < 1) {
        logger.severe(
            format(
                " FFX shutting down: end range for -fR must be at least 1; input range %d to %d",
                forceResiduesStart, forceResiduesEnd));
      }
    }

    if (group.algorithm != 5) {
      if (!residueGroup.listResidues.equalsIgnoreCase("none")) {
        StringBuilder info = new StringBuilder("\n Evaluating rotamers for residues ");
        for (String i : resList) {
          info.append(format("%s, ", i));
        }
        logger.info(info.toString());
      } else if (allStartResID == -1) {
        logger.info("\n Evaluating rotamers for residues " + residueGroup.start + " to "
            + residueGroup.finish);
      } else {
        logger.info("\n Evaluating rotamers for all residues beginning at " + allStartResID);
      }
    } else {
      if (!residueGroup.listResidues.equalsIgnoreCase("none")) {
        StringBuilder info = new StringBuilder("\n Evaluating rotamers for boxes with residues ");
        for (String i : resList) {
          info.append(format("%s, ", i));
        }
        logger.info(info.toString());
      } else if (allStartResID == -1) {
        logger.info("\n Evaluating rotamers for boxes " + (boxStart + 1) + " to " + (boxEnd + 1));
      } else {
        logger.info("\n Evaluating rotamers for all boxes beginning at " + (boxStart + 1));
      }
    }
  }

  /**
   * addListResidues.
   *
   * @param resList a {@link java.util.List} object.
   */
  private void addListResidues(List<String> resList) {
    if (!residueGroup.listResidues.equalsIgnoreCase("none")) {
      String[] tok = residueGroup.listResidues.split(",");
      for (String t : tok) {
        logger.info(" Adding " + t);
        resList.add(t);
      }
    }
  }

  /** Sets the standard values for properties in rotamer optimization. */
  private void setRotOptProperties() {
    // General
    rotamerOptimization.setTwoBodyCutoff(energyGroup.twoBodyCutoff);
    rotamerOptimization.setThreeBodyCutoff(energyGroup.threeBodyCutoff);
    rotamerOptimization.setThreeBodyEnergy(energyGroup.threeBody);
    rotamerOptimization.setUseGoldstein(!group.dee);
    rotamerOptimization.setRevert(group.revert);
    rotamerOptimization.setPruning(energyGroup.prune);
    rotamerOptimization.setDistanceCutoff(energyGroup.cutoff);
    boolean monteCarloBool = false;
    if (group.monteCarlo > 1) {
      monteCarloBool = true;
    }
    rotamerOptimization.setMonteCarlo(monteCarloBool, group.monteCarlo);

    File energyRestartFile = null;
    if (!group.energyRestart.equalsIgnoreCase("none")) {
      energyRestartFile = new File(group.energyRestart);
      rotamerOptimization.setEnergyRestartFile(energyRestartFile);
    }

    // Window
    if (group.algorithm == 4) {
      rotamerOptimization.setWindowSize(windowGroup.window);
      rotamerOptimization.setIncrement(windowGroup.increment);
      rotamerOptimization.setForcedResidues(forceResiduesStart, forceResiduesEnd);
    }

    // Box
    if (group.algorithm == 5) {
      if (boxGroup.approxBoxLength < 0) {
        logger.info(" Negative box length value changed to -1 * input.");
        boxGroup.approxBoxLength *= -1;
      }
      rotamerOptimization.setBoxBorderSize(boxGroup.boxBorderSize);
      rotamerOptimization.setApproxBoxLength(boxGroup.approxBoxLength);
      rotamerOptimization.setNumXYZBoxes(numXYZBoxes);
      rotamerOptimization.setBoxInclusionCriterion(boxGroup.boxInclusionCriterion);
    }
  }

  /**
   * Choices are independent residues (1), all with rotamer elimination (2), all brute force (3),
   * sliding window (4), or box optimization (5).
   *
   * @return Returns the algorithm choice.
   */
  public int getAlgorithm() {
    return group.algorithm;
  }

  public void setAlgorithm(int algorithm) {
    group.algorithm = algorithm;
  }

  /**
   * Ponder and Richards (1) or Richardson (2) rotamer library.
   *
   * @return Returns the Rotamer library.
   */
  public int getLibrary() {
    return group.library;
  }

  public void setLibrary(int library) {
    group.library = library;
  }

  /**
   * Nucleic acid library: currently only Richardson available.
   *
   * @return Returns a String for the nucleic acid library.
   */
  public String getNaLibraryName() {
    return group.naLibraryName;
  }

  public void setNaLibraryName(String naLibraryName) {
    group.naLibraryName = naLibraryName;
  }

  /**
   * Use dead-end elimination criteria instead of Goldstein criteria.
   *
   * @return Returns true if using DEE instead of Goldstein.
   */
  public boolean isDee() {
    return group.dee;
  }

  public void setDee(boolean dee) {
    group.dee = dee;
  }

  /**
   * Single character chain ID of the residues to optimize.
   *
   * @return Returns the Chain name.
   */
  public String getChain() {
    return residueGroup.chain;
  }

  public void setChain(String chain) {
    residueGroup.chain = chain;
  }

  /**
   * Starting residue to perform the optimization on (-1 exits). For box optimization, first box to
   * optimize.
   *
   * @return Returns the starting index.
   */
  public int getStart() {
    return residueGroup.start;
  }

  public void setStart(int start) {
    residueGroup.start = start;
  }

  /**
   * Final residue to perform the optimization on (-1 exits). For box optimization, final box to
   * optimize.
   *
   * @return Returns the finish index.
   */
  public int getFinish() {
    return residueGroup.finish;
  }

  public void setFinish(int finish) {
    residueGroup.finish = finish;
  }

  /**
   * Cutoff distance for two-body interactions.
   *
   * @return Returns the 2-body cutoff.
   */
  public double getTwoBodyCutoff() {
    return energyGroup.twoBodyCutoff;
  }

  public void setTwoBodyCutoff(double twoBodyCutoff) {
    energyGroup.twoBodyCutoff = twoBodyCutoff;
  }

  /**
   * -T or --threeBody Include 3-Body interactions in the elimination criteria.
   *
   * @return Returns true if 3-body interactions are being used.
   */
  public boolean isThreeBody() {
    return energyGroup.threeBody;
  }

  public void setThreeBody(boolean threeBody) {
    energyGroup.threeBody = threeBody;
  }

  /**
   * Cutoff distance for three-body interactions.
   *
   * @return Returns the 3-body cutoff.
   */
  public double getThreeBodyCutoff() {
    return energyGroup.threeBodyCutoff;
  }

  public void setThreeBodyCutoff(double threeBodyCutoff) {
    energyGroup.threeBodyCutoff = threeBodyCutoff;
  }

  /**
   * Prune no clashes (0), only single clashes (1), or all clashes (2).
   *
   * @return Returns the pruning condition.
   */
  public int getPrune() {
    return energyGroup.prune;
  }

  public void setPrune(int prune) {
    energyGroup.prune = prune;
  }

  /**
   * Optimize all residues beginning from the passed value (overrides other options). for box
   * optimization, optimizes all boxes beginning from the passed index. Default is to optimize all
   * residues.
   *
   * @return Returns the residue / box to start from.
   */
  public int getAll() {
    return residueGroup.all;
  }

  public void setAll(int all) {
    residueGroup.all = all;
  }

  /**
   * Revert unfavorable changes.
   *
   * @return Returns true if unfavorable changes are reverted.
   */
  public boolean isRevert() {
    return group.revert;
  }

  public void setRevert(boolean revert) {
    group.revert = revert;
  }

  /**
   * Energy restart file from a previous run (requires that all parameters are the same).
   *
   * @return Returns the energy restart file to use.
   */
  public String getEnergyRestart() {
    return group.energyRestart;
  }

  public void setEnergyRestart(String energyRestart) {
    group.energyRestart = energyRestart;
  }

  /**
   * Do not include starting coordinates as their own rotamer.
   *
   * @return Returns true if original side-chain coordinates should not be used as a rotamer.
   */
  public boolean isNoOriginal() {
    return group.noOriginal;
  }

  public void setNoOriginal(boolean noOriginal) {
    group.noOriginal = noOriginal;
  }

  /**
   * -E or --decompose Print energy decomposition for the input structure (no optimization).
   *
   * @return Returns true if the input structure should undergo an energy decomposition.
   */
  public boolean isDecompose() {
    return group.decompose;
  }

  public void setDecompose(boolean decompose) {
    group.decompose = decompose;
  }

  /**
   * Choose a list of individual residues to optimize (eg. A11,A24,B40).
   *
   * @return Returns the list of selected residues.
   */
  public String getListResidues() {
    return residueGroup.listResidues;
  }

  public void setListResidues(String listResidues) {
    residueGroup.listResidues = listResidues;
  }

  /**
   * Follow elimination criteria with 'n' Monte Carlo steps, or enumerate all remaining
   * conformations, whichever is smaller.
   *
   * @return Returns the number of Monte Carlo optimization steps to apply.
   */
  public int getMonteCarlo() {
    return group.monteCarlo;
  }

  public void setMonteCarlo(int monteCarlo) {
    group.monteCarlo = monteCarlo;
  }

  /**
   * Save eliminated singles and eliminated pairs to a text file (global and box optimization).
   *
   * @return Returns true to Save eliminated rotamers to a file.
   */
  public boolean isSaveOutput() {
    return group.saveOutput;
  }

  public void setSaveOutput(boolean saveOutput) {
    group.saveOutput = saveOutput;
  }

  /**
   * Size of the sliding window with respect to adjacent residues (default = 7).
   *
   * @return Returns the sliding window size.
   */
  public int getWindow() {
    return windowGroup.window;
  }

  public void setWindow(int window) {
    windowGroup.window = window;
  }

  /**
   * Sliding window increment (default = 3).
   *
   * @return Returns the sliding window increment.
   */
  public int getIncrement() {
    return windowGroup.increment;
  }

  public void setIncrement(int increment) {
    windowGroup.increment = increment;
  }

  /**
   * The sliding window and box cutoff radius (Angstroms).
   *
   * @return Returns the sliding window cutoff radius.
   */
  public double getCutoff() {
    return energyGroup.cutoff;
  }

  public void setCutoff(double cutoff) {
    energyGroup.cutoff = cutoff;
  }

  /**
   * The threshold for pruning clashes. If two self-energies on the same residue have an energy
   * difference greater than 25 kcal/mol, the high energy rotamers get pruned.
   *
   * @return Returns the clash threshold for self energies.
   */
  public double getClashThreshold() {
    return energyGroup.clashThreshold;
  }

  public void setClashThreshold(double clashThreshold) {
    energyGroup.clashThreshold = clashThreshold;
  }

  /**
   * The threshold for pruning clashes. If two pair-energies on the same residues have an energy
   * difference greater than 25 kcal/mol, the high energy rotamers get pruned.
   *
   * @return Returns the clash threshold for pair energies.
   */
  public double getPairClashThreshold() {
    return energyGroup.pairClashThreshold;
  }

  public void setPairClashThreshold(double pairClashThreshold) {
    energyGroup.pairClashThreshold = pairClashThreshold;
  }

  /**
   * Force residues in this range to be considered for sliding window radii, regardless of whether
   * they lack rotamers.
   *
   * @return Returns forced residues.
   */
  public String getForceResidues() {
    return group.forceResidues;
  }

  public void setForceResidues(String forceResidues) {
    group.forceResidues = forceResidues;
  }

  /**
   * The number of boxes along X, Y, and Z (default: '3,3,3').
   *
   * @return Returns the number of boxes.
   */
  public String getNumBoxes() {
    return boxGroup.numBoxes;
  }

  public void setNumBoxes(String numBoxes) {
    boxGroup.numBoxes = numBoxes;
  }

  /**
   * Extent of overlap between optimization boxes (default: 0.0 A).
   *
   * @return Returns the overlap between optimization boxes.
   */
  public double getBoxBorderSize() {
    return boxGroup.boxBorderSize;
  }

  public void setBoxBorderSize(double boxBorderSize) {
    boxGroup.boxBorderSize = boxBorderSize;
  }

  /**
   * Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes). Box sizes are
   * rounded up to make a whole number of boxes along each axis (default of 0 disables this
   * function).
   *
   * @return Returns the approximate box length.
   */
  public double getApproxBoxLength() {
    return boxGroup.approxBoxLength;
  }

  public void setApproxBoxLength(double approxBoxLength) {
    boxGroup.approxBoxLength = approxBoxLength;
  }

  /**
   * Criterion to use for adding residues to boxes. (1) uses C alpha only (N1/9 for nucleic acids)
   * (2) uses any atom. (3) uses any rotamer.
   *
   * @return Returns the Box inclusion criteria.
   */
  public int getBoxInclusionCriterion() {
    return boxGroup.boxInclusionCriterion;
  }

  public void setBoxInclusionCriterion(int boxInclusionCriterion) {
    boxGroup.boxInclusionCriterion = boxInclusionCriterion;
  }

  /**
   * Collection of ManyBody Options.
   */
  private static class ManyBodyOptionGroup {

    /**
     * -a or --algorithm Choices are independent residues (1), all with rotamer elimination (2), all
     * brute force (3), sliding window (4), or box optimization (5).
     */
    @Option(
        names = {"-a", "--algorithm"},
        paramLabel = "0",
        defaultValue = "0",
        description =
            "Algorithm: default automatic settings (0), independent residues (1), all with rotamer elimination (2), all brute force (3), sliding window (4), or box optimization (5)")
    private int algorithm;

    /** -L or --library Choose either Ponder and Richards (1) or Richardson (2) rotamer library. */
    @Option(
        names = {"-L", "--library"},
        paramLabel = "2",
        defaultValue = "2",
        description = "Ponder and Richards (1) or Richardson (2) rotamer library.")
    private int library;

    /** -Ln or --libraryNucleic Choose a nucleic acid library: currently only Richardson available. */
    @Option(
        names = {"--Ln", "--libraryNucleic"},
        paramLabel = "Richardson",
        defaultValue = "Richardson",
        description = "Nucleic acid library to select: [Richardson]")
    private String naLibraryName;

    /** --dee or --deadEnd Use dead-end elimination criteria instead of Goldstein criteria. */
    @Option(
        names = {"--dee", "--deadEnd"},
        paramLabel = "false",
        defaultValue = "false",
        description = "Use dead-end elimination criteria instead of Goldstein criteria.")
    private boolean dee;

    /** -z or --revert Revert unfavorable changes. */
    @Option(
        names = {"-z", "--revert"},
        defaultValue = "true",
        description = "Revert unfavorable changes.")
    private boolean revert;

    /**
     * --eR or --energyRestart Load energy restart file from a previous run (requires that all
     * parameters are the same).
     */
    @Option(
        names = {"--eR", "--energyRestart"},
        paramLabel = "none",
        defaultValue = "none",
        description =
            "Load energy restart file from a previous run (requires that all parameters are the same).")
    private String energyRestart;

    /** -o or --noOriginal Do not include starting coordinates as their own rotamer. */
    @Option(
        names = {"-O", "--noOriginal"},
        defaultValue = "false",
        description = "Do not include starting coordinates as their own rotamer.")
    private boolean noOriginal;

    /** -E or --decompose Print energy decomposition for the input structure (no optimization). */
    @Option(
        names = {"-E", "--decompose"},
        defaultValue = "false",
        description = "Print energy decomposition for the input structure (no optimization!).")
    private boolean decompose;

    /**
     * --pH or --titrationPH Optimize the titration state of ASP, GLU, HIS and LYS residues.
     */
    @Option(
        names = {"--pH", "--titrationPH"},
        paramLabel = "0",
        defaultValue = "0",
        description = " Optimize the titration state of ASP, GLU, HIS and LYS residues at the given pH (pH = 0 turns off titration")
    private double titrationPH;

    /**
     * --mC or --monteCarlo Follow elimination criteria with 'n' Monte Carlo steps, or enumerate all
     * remaining conformations, whichever is smaller.
     */
    @Option(
        names = {"--mC", "--monteCarlo"},
        paramLabel = "-1",
        defaultValue = "-1",
        description =
            "Follow elimination criteria with (n) Monte Carlo steps, or enumerate all remaining conformations, whichever is smaller.")
    private int monteCarlo;

    /**
     * -out or --output Save eliminated singles and eliminated pairs to a text file (global and box
     * optimization).
     */
    //  @Option(
    //      names = {"--out", "--output"},
    //      defaultValue = "false",
    //      description = "Save eliminated singles and eliminated pairs to a text file.")
    private boolean saveOutput;

    /**
     * -fR or --forceResidues Force residues in this range to be considered for sliding window radii,
     * regardless of whether they lack rotamers.
     */
    @Option(
        names = {"--fR", "--forceResidues"},
        paramLabel = "-1,-1",
        defaultValue = "-1,-1",
        description =
            "Force residues in this range to be considered for sliding window radii, regardless of whether they lack rotamers.")
    private String forceResidues;

  }

  /**
   * Collection of ManyBody Box Optimization Options.
   */
  private static class ManyBodyBoxOptionGroup {

    /** -nB or --numBoxes Specify number of boxes along X, Y, and Z (default: '3,3,3'). */
    @Option(
        names = {"--nB", "--numBoxes"},
        paramLabel = "3,3,3",
        defaultValue = "3,3,3",
        description = "Specify number of boxes along X, Y, and Z (default: 3,3,3)")
    private String numBoxes;

    /** -bB or --boxBorderSize Extent of overlap between optimization boxes (default: 0.0 A). */
    @Option(
        names = {"--bB", "--boxBorderSize"},
        paramLabel = "0.0",
        defaultValue = "0.0",
        description = "Extent of overlap between optimization boxes in Angstroms.")
    private double boxBorderSize;

    /**
     * -bL or --approxBoxLength Approximate side lengths of boxes to be constructed (over-rides
     * numXYZBoxes). Box sizes are rounded up to make a whole number of boxes along each axis
     * (default of 0 disables this function).
     */
    @Option(
        names = {"--bL", "--approxBoxLength"},
        paramLabel = "20.0",
        defaultValue = "20.0",
        description = "Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes).")
    private double approxBoxLength;

    /**
     * -bC or --boxInclusionCriterion Criterion to use for adding residues to boxes. (1) uses C alpha
     * only (N1/9 for nucleic acids) (2) uses any atom. (3) uses any rotamer
     */
    @Option(
        names = {"--bC", "--boxInclusionCriterion"},
        paramLabel = "1",
        defaultValue = "1",
        description =
            "Criterion to use for adding a residue to a box: (1) uses C alpha only (N1/9 for nucleic acids), (2) uses any atom, and (3) uses any rotamer")
    private int boxInclusionCriterion;
  }

  /**
   * Collection of ManyBody Window Optimization Options.
   */
  private static class ManyBodyWindowOptionGroup {

    /** --window Size of the sliding window with respect to adjacent residues (default = 7). */
    @Option(
        names = {"--window"},
        paramLabel = "7",
        defaultValue = "7",
        description = "Size of the sliding window with respect to adjacent residues.")
    private int window;

    /** --increment Sliding window increment (default = 3). */
    @Option(
        names = {"--increment"},
        paramLabel = "3",
        defaultValue = "3",
        description = "Sliding window increment.")
    private int increment;

  }

  /**
   * Collection of ManyBody Energy Optimization Options.
   */
  private static class ManyBodyEnergyOptionGroup {

    /** --radius The sliding window and box cutoff radius (Angstroms). */
    @Option(
        names = {"--radius"},
        paramLabel = "2.0",
        defaultValue = "2.0",
        description = "The sliding box and window cutoff radius (Angstroms).")
    private double cutoff;

    /** --tC or --twoBodyCutoff Cutoff distance for two-body interactions. */
    @Option(
        names = {"--tC", "--twoBodyCutoff"},
        paramLabel = "3.0",
        defaultValue = "3.0",
        description = "Cutoff distance for two body interactions.")
    private double twoBodyCutoff;

    /** -T or --threeBody Include 3-Body interactions in the elimination criteria. */
    @Option(
        names = {"-T", "--threeBody"},
        defaultValue = "false",
        description = "Include 3-Body interactions in the elimination criteria.")
    private boolean threeBody;

    /** --thC or --threeBodyCutoff Cutoff distance for three-body interactions. */
    @Option(
        names = {"--thC", "--threeBodyCutoff"},
        paramLabel = "3.0",
        defaultValue = "3.0",
        description = "Cutoff distance for three-body interactions.")
    private double threeBodyCutoff;

    /** --pr or --prune Prune no clashes (0), only single clashes (1), or all clashes (2). */
    @Option(
        names = {"--pr", "--prune"},
        paramLabel = "1",
        defaultValue = "1",
        description = "Prune no clashes (0), only single clashes (1), or all clashes (2)")
    private int prune;

    /**
     * --clashThreshold The threshold for pruning clashes. If two self-energies on the same residue
     * have an energy difference greater than 25 kcal/mol, the high energy rotamers get pruned.
     */
    @Option(
        names = {"--clashThreshold"},
        paramLabel = "25.0",
        defaultValue = "25.0",
        description = "The threshold for pruning clashes.")
    private double clashThreshold;

    /**
     * --pairClashThreshold The threshold for pruning clashes. If two pair-energies on the same
     * residues have an energy difference greater than 25 kcal/mol, the high energy rotamers get
     * pruned.
     */
    @Option(
        names = {"--pairClashThreshold"},
        paramLabel = "25.0",
        defaultValue = "25.0",
        description = "The threshold for pruning pair clashes.")
    private double pairClashThreshold;
  }

  /**
   * Collection of ManyBody Residue Selection Options.
   */
  private static class ManyBodyResidueOptionGroup {

    /** --ch or --chain Single character chain ID of the residues to optimize. */
    @Option(
        names = {"--ch", "--chain"},
        paramLabel = "-1",
        defaultValue = "-1",
        description = "Single character chain ID of the residues to optimize.")
    private String chain;

    /**
     * -s or --start Starting residue to perform the optimization on (-1 exits). For box
     * optimization, first box to optimize.
     */
    @Option(
        names = {"-s", "--start"},
        paramLabel = "-1",
        defaultValue = "-1",
        description =
            "Starting residue to perform the optimization on (-1 exits). For box optimization, first box to optimize.")
    private int start;

    /**
     * --fi or --final Final residue to perform the optimization on (-1 exits). For box optimization,
     * final box to optimize.
     */
    @Option(
        names = {"--fi", "--final"},
        paramLabel = "-1",
        defaultValue = "-1",
        description =
            "Final residue to perform the optimization on (-1 exits). For box optimization, final box to optimize.")
    private int finish;

    /**
     * -x or --all Optimize all residues beginning from the passed value (overrides other options);
     * for box optimization, optimizes all boxes beginning from the passed index. Default is to
     * optimize all residues.
     */
    @Option(
        names = {"-x", "--all"},
        paramLabel = "-1",
        defaultValue = "-1",
        description =
            "Optimize all residues beginning from the passed value (overrides other options); for box optimization, optimizes all boxes beginning from the passed index.")
    private int all;

    /** --lR or --listResidues Choose a list of individual residues to optimize (eg. A11,A24,B40). */
    @Option(
        names = {"--lR", "--listResidues"},
        paramLabel = "none",
        defaultValue = "none",
        description = "Choose a list of individual residues to optimize (eg. A11,A24,B40).")
    private String listResidues;

  }

  //    /**
  //     * Saves all eliminated rotamers to an ouput file called "eliminated.csv"
  //     * when the many body command is run with the following syntax and flags:
  //     * ffxc ManyBody --out ... file.pdb &gt;&gt; file.log.
  //     *
  //     * @throws java.io.IOException Throws an exception when output is non piped to a log
  //     *                             file. The --out flag relies on the presence of a log file
  // where output is
  //     *                             piped.
  //     */
  //    public void saveEliminatedRotamers() throws IOException {
  //        if (saveOutput) {
  //            rotamerOptimization.outputEliminated();
  //        }
  //    }

}
