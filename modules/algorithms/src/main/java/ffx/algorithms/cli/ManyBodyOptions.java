// ******************************************************************************
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
// ******************************************************************************
package ffx.algorithms.cli;

import static java.lang.Integer.parseInt;

import ffx.algorithms.optimize.RotamerOptimization;
import ffx.algorithms.optimize.RotamerOptimization.Algorithm;
import ffx.numerics.math.DoubleMath;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.*;

import java.io.File;
import java.util.*;
import java.util.logging.Logger;

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
public class ManyBodyOptions {

  private static final Logger logger = Logger.getLogger(ManyBodyOptions.class.getName());

  /**
   * The ArgGroup keeps the ManyBodyOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Many-Body Optimization Options%n", validate = false)
  private final ManyBodyOptionGroup group = new ManyBodyOptionGroup();

  /**
   * The ArgGroup keeps the BoxOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Many-Body Box Optimization Options%n", validate = false)
  private final BoxOptionGroup boxGroup = new BoxOptionGroup();

  /**
   * The ArgGroup keeps the WindowOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Many-Body Window Optimization Options%n", validate = false)
  private final WindowOptionGroup windowGroup = new WindowOptionGroup();

  /**
   * The ArgGroup keeps the WindowOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Many-Body Energy Expansion and Cut-off Options%n", validate = false)
  private final EnergyOptionGroup energyGroup = new EnergyOptionGroup();

  /**
   * The ArgGroup keeps the ResidueOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n Many-Body Residue Selection Options%n", validate = false)
  private final ResidueOptionGroup residueGroup = new ResidueOptionGroup();

  private RotamerOptimization rotamerOptimization;
  private RotamerLibrary rotamerLibrary;

  /**
   * initRotamerOptimization.
   *
   * @param rotamerOptimization a {@link RotamerOptimization} object.
   * @param activeAssembly a {@link ffx.potential.MolecularAssembly} object.
   */
  public void initRotamerOptimization(RotamerOptimization rotamerOptimization,
      MolecularAssembly activeAssembly) {
    this.rotamerOptimization = rotamerOptimization;

    // Collect the residues to optimize.
    List<Residue> residues = collectResidues(activeAssembly);
    rotamerOptimization.setResidues(residues);
    rotamerOptimization.setRotamerLibrary(rotamerLibrary);

    // If the user has not selected an algorithm, it will be chosen based on the number of residues.
    Algorithm algorithm = getAlgorithm(residues.size());

    // Configure general options.
    rotamerOptimization.setDecomposeOriginal(group.decompose);
    rotamerOptimization.setUseGoldstein(!group.dee);
    rotamerOptimization.setRevert(group.revert);
    boolean monteCarloBool = group.monteCarlo > 1;
    rotamerOptimization.setMonteCarlo(monteCarloBool, group.monteCarlo);
    File energyRestartFile;
    if (!group.energyRestart.equalsIgnoreCase("none")) {
      energyRestartFile = new File(group.energyRestart);
      rotamerOptimization.setEnergyRestartFile(energyRestartFile);
    }

    // Configure Energy Expansion and Pruning Options.
    rotamerOptimization.setTwoBodyCutoff(energyGroup.twoBodyCutoff);
    rotamerOptimization.setThreeBodyEnergy(energyGroup.threeBody);
    rotamerOptimization.setThreeBodyCutoff(energyGroup.threeBodyCutoff);
    rotamerOptimization.setDistanceCutoff(energyGroup.cutoff);
    rotamerOptimization.setPruning(energyGroup.prune);
    rotamerOptimization.setSingletonClashThreshold(energyGroup.clashThreshold);
    rotamerOptimization.setPairClashThreshold(energyGroup.pairClashThreshold);

    // Window
    if (algorithm == Algorithm.WINDOW) {
      rotamerOptimization.setWindowSize(windowGroup.window);
      rotamerOptimization.setIncrement(windowGroup.increment);
    } else if (algorithm == Algorithm.BOX) {
      // Box
      parseBoxSelection();
      if (boxGroup.approxBoxLength < 0) {
        logger.info(" Negative box length value changed to -1 * input.");
        boxGroup.approxBoxLength *= -1;
      }
      rotamerOptimization.setBoxBorderSize(boxGroup.boxBorderSize);
      rotamerOptimization.setApproxBoxLength(boxGroup.approxBoxLength);
      rotamerOptimization.setNumXYZBoxes(boxGroup.numXYZBoxes);
      rotamerOptimization.setBoxInclusionCriterion(boxGroup.boxInclusionCriterion);
      rotamerOptimization.setBoxStart(boxGroup.initialBox);
      rotamerOptimization.setBoxEnd(boxGroup.finalBox);
      rotamerOptimization.setTitrationBoxes(boxGroup.boxTitration);
      rotamerOptimization.setTitrationBoxSize(boxGroup.approxBoxLength);
    }
  }

  /**
   * Collect residues based on residue selection flags.
   *
   * @param activeAssembly a {@link ffx.potential.MolecularAssembly} object.
   */
  public List<Residue> collectResidues(MolecularAssembly activeAssembly) {

    // Force re-initialization RotamerLibrary prior to collecting residues and rotamers.
    initRotamerLibrary(true);

    // First, interpret the residueGroup.listResidues flag if its set.
    if (!residueGroup.listResidues.isEmpty() && !residueGroup.listResidues.equalsIgnoreCase("none")) {
      List<String> stringList = new ArrayList<>();
      String[] tok = residueGroup.listResidues.split(",");
      Collections.addAll(stringList, tok);
      List<Residue> residueList = new ArrayList<>();
      Polymer[] polymers = activeAssembly.getChains();
      for (String s : stringList) {
        Character chainID = s.charAt(0);
        int i = parseInt(s.substring(1));
        for (Polymer polymer : polymers) {
          if (polymer.getChainID() == chainID) {
            List<Residue> residues = polymer.getResidues();
            for (Residue residue : residues) {
              if (residue.getResidueNumber() == i) {
                Rotamer[] rotamers = residue.setRotamers(rotamerLibrary);
                if (rotamers != null && rotamers.length > 0) {
                  residueList.add(residue);
                }
              }
            }
          }
        }
      }
      return residueList;
    }

    // Check that the finish flag is greater than the start flag.
    if (residueGroup.finish < residueGroup.start) {
      residueGroup.finish = Integer.MAX_VALUE;
    }
    Character chainID = null;
    if (!residueGroup.chain.equalsIgnoreCase("-1")) {
      chainID = residueGroup.chain.charAt(0);
    }

    // Otherwise, collect all residues with a rotamer.
    List<Residue> residueList = new ArrayList<>();
    Polymer[] polymers = activeAssembly.getChains();
    for (Polymer polymer : polymers) {
      // Enforce requested chainID.
      if (chainID != null && chainID != polymer.getChainID()) {
        continue;
      }
      List<Residue> residues = polymer.getResidues();
      for (Residue residue : residues) {
        int resID = residue.getResidueNumber();
        // Enforce requested residue range.
        if (resID >= residueGroup.start && resID <= residueGroup.finish) {
          Rotamer[] rotamers = residue.setRotamers(rotamerLibrary);
          if (rotamers != null) {
            residueList.add(residue);
          }
        }
      }
    }

    return residueList;
  }

  /**
   * Returns the user selected algorithm or one chosen based on number of residues.
   *
   * @return Returns the algorithm choice.
   */
  public Algorithm getAlgorithm(int numResidues) {
    if (group.algorithm == 0) {
      if (numResidues < 100) {
        return Algorithm.ALL;
      } else {
        return Algorithm.BOX;
      }
    }
    return Algorithm.getAlgorithm(group.algorithm);
  }

  public boolean getUsingOriginalCoordinates() {
    return !group.noOriginal;
  }

  public void setOriginalCoordinates(boolean useOrig) {
    group.noOriginal = !useOrig;
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

  public RotamerLibrary getRotamerLibrary(boolean reinit) {
    initRotamerLibrary(reinit);
    return rotamerLibrary;
  }

  private void initRotamerLibrary(boolean reinit) {
    if (rotamerLibrary == null || reinit) {
      boolean useOrigCoordsRotamer = !group.noOriginal;
      if (group.decompose) {
        useOrigCoordsRotamer = true;
      }
      rotamerLibrary = new RotamerLibrary(
          RotamerLibrary.ProteinLibrary.intToProteinLibrary(group.library), useOrigCoordsRotamer);
    }
  }


  /** Set allStartResID, boxStart and boxEnd */
  private void parseBoxSelection() {
    // Parse the numBoxes flag.
    String input = boxGroup.numBoxes;
    Scanner boxNumInput = new java.util.Scanner(input);
    boxNumInput.useDelimiter(",");
    int inputLoopCounter = 0;
    int[] numXYZBoxes = new int[3];
    numXYZBoxes[0] = 3; // Default
    while (inputLoopCounter < 3) {
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

    // Initialize dimensions not provided.
    for (int i = inputLoopCounter; i < 3; i++) {
      numXYZBoxes[i] = numXYZBoxes[0];
    }

    // Correct input errors.
    int totalCount = 1;
    for (int i = 0; i < 3; i++) {
      if (numXYZBoxes[i] == 0) {
        numXYZBoxes[i] = 3;
        logger.info(" Input of 0 to nB reset to default of 3.");
      } else if (numXYZBoxes[i] < 0) {
        numXYZBoxes[i] = -1 * numXYZBoxes[i];
        logger.info(" Input of negative number to nB reset to positive number");
      }
      totalCount *= numXYZBoxes[i];
    }

    boxGroup.numXYZBoxes = numXYZBoxes;

    if (boxGroup.initialBox < 0) {
      boxGroup.initialBox = 0;
    }
    if (boxGroup.finalBox < 0 || boxGroup.finalBox > totalCount) {
      boxGroup.finalBox = totalCount;
    }

  }

  /** Sets the standard values for properties in rotamer optimization. */
  private void setRotOptProperties(Algorithm algorithm) {

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

  public boolean getOnlyTitration() {
    return residueGroup.onlyTitration;
  }

  public void setOnlyTitration(boolean onlyTitration) {
    residueGroup.onlyTitration = onlyTitration;
  }

  public int getInterestedResidue() {
    return residueGroup.interestedResidue;
  }

  public void setInterestedResidue(int interestedResidue) {
    residueGroup.interestedResidue = interestedResidue;
  }

  public double getInclusionCutoff() {
    return residueGroup.inclusionCutoff;
  }

  public void setInclusionCutoff(double inclusionCutoff) {
    residueGroup.inclusionCutoff = inclusionCutoff;
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

  public void setBoxTitration(boolean boxTitration){boxGroup.boxTitration = boxTitration;}

  public boolean getBoxTitration(){return boxGroup.boxTitration;}

  public void setTitrationPH(double pH) {
    group.titrationPH = pH;
  }

  public double getTitrationPH() {
    return group.titrationPH;
  }

  public void setPHRestraint(double pHRestraint) {
    energyGroup.pHRestraint = pHRestraint;
  }

  public double getPHRestraint() {
    return energyGroup.pHRestraint;
  }

  public boolean isTitrating() {
    return group.titrationPH == 0;
  }

  public boolean getTitration() {
    return group.titration;
  }

  public String selectInclusionResidues(final List<Residue> residueList, int mutatingResidue, boolean onlyTitration,
                                       double inclusionCutoff){
    String listResidues = "";
    if (mutatingResidue != -1 && inclusionCutoff != -1) {
      List<Integer> residueNumber = new ArrayList<>();
      for (Residue residue : residueList) {
        residueNumber.add(residue.getResidueNumber());
      }
      double[] mutatingResCoor = new double[3];
      int index = residueNumber.indexOf(mutatingResidue);
      mutatingResCoor = residueList.get(index).getAtomByName("CA", true).getXYZ(mutatingResCoor);
      for (Residue residue: residueList) {
        double[] currentResCoor = new double[3];
        currentResCoor = residue.getAtomByName("CA", true).getXYZ(currentResCoor);
        double dist = DoubleMath.dist(mutatingResCoor, currentResCoor);
        if (dist < inclusionCutoff) {
          listResidues += "," + residue.getChainID() + residue.getResidueNumber();
        }
      }
      listResidues = listResidues.substring(1);
    } else if (onlyTitration){
      String[] titratableResidues = new String[]{"HIS", "HIE", "HID", "GLU", "GLH", "ASP", "ASH", "LYS", "LYD", "CYS", "CYD"};
      List<String> titratableResiudesList = Arrays.asList(titratableResidues);
      for (Residue residue : residueList) {
        if (titratableResiudesList.contains(residue.getName())) {
          String titrateResNum = Integer.toString(residue.getResidueNumber());
          if(!listResidues.contains(titrateResNum)){
            listResidues += "," + residue.getChainID()+ titrateResNum;
          }
          if (inclusionCutoff != -1){
            for (Residue residue2: residueList) {
              boolean includeResidue = evaluateAllRotDist(residue, residue2, inclusionCutoff);
              if(includeResidue){
                String residue2Number = Integer.toString(residue2.getResidueNumber());
                if(!listResidues.contains(residue2Number)){
                  listResidues += "," + residue2.getChainID()+ residue2Number;
                }
              }
            }
          }

        }

      }

      listResidues = listResidues.substring(1);
    }
    return listResidues;
  }
  private static boolean evaluateAllRotDist(Residue residueA, Residue residueB, double inclusionCutoff){
    residueA.setRotamers(RotamerLibrary.getDefaultLibrary());
    residueB.setRotamers(RotamerLibrary.getDefaultLibrary());
    Rotamer[] rotamersA = residueA.getRotamers();
    Rotamer[] rotamersB = residueB.getRotamers();
    double[] aCoor = new double[3];
    double[] bCoor = new double[3];
    try {
      int a = rotamersA.length;
      int b = rotamersB.length;
    } catch (Exception e){
      return false;
    }

    for(Rotamer rotamerA: rotamersA){
      residueA.setRotamer(rotamerA);
      for(Rotamer rotamerB: rotamersB){
        residueB.setRotamer(rotamerB);
        for(Atom atomA: residueA.getAtomList()){
          for(Atom atomB: residueB.getAtomList()){
            double dist = DoubleMath.dist(atomA.getXYZ(aCoor), atomB.getXYZ(bCoor));
            if(dist <= inclusionCutoff){
              return true;
            }
          }
        }
      }
    }

    return false;
  }


  /**
   * Collection of ManyBody Options.
   */
  private static class ManyBodyOptionGroup {

    /**
     * -a or --algorithm Choices are independent residues (1), all with rotamer elimination (2), all
     * brute force (3), sliding window (4), or box optimization (5).
     */
    @Option(names = {"-a",
        "--algorithm"}, paramLabel = "0", defaultValue = "0", description = "Algorithm: default automatic settings (0), independent residues (1), all with rotamer elimination (2), all brute force (3), sliding window (4), or box optimization (5)")
    private int algorithm;

    /** -L or --library Choose either Ponder and Richards (1) or Richardson (2) rotamer library. */
    @Option(names = {"-L",
        "--library"}, paramLabel = "2", defaultValue = "2", description = "Ponder and Richards (1) or Richardson (2) rotamer library.")
    private int library;

    /** -nl or --nucleicLibrary Choose a nucleic acid library: currently only Richardson available. */
    // @Option(
    //    names = {"--nl", "--nucleiclibrary"},
    //    paramLabel = "Richardson",
    //    defaultValue = "Richardson",
    //    description = "Nucleic acid library to select: [Richardson]")
    private String naLibraryName = "Richardson";

    /** --dee or --deadEnd Use dead-end elimination criteria instead of Goldstein criteria. */
    @Option(names = {"--dee",
        "--deadEnd"}, paramLabel = "false", defaultValue = "false", description = "Use dead-end elimination criteria instead of Goldstein criteria.")
    private boolean dee;

    /** -z or --revert Revert unfavorable changes. */
    @Option(names = {"-z",
        "--revert"}, defaultValue = "true", description = "Revert unfavorable changes.")
    private boolean revert;

    /**
     * --eR or --energyRestart Load energy restart file from a previous run (requires that all
     * parameters are the same).
     */
    @Option(names = {"--eR",
        "--energyRestart"}, paramLabel = "none", defaultValue = "none", description = "Load energy restart file from a previous run (requires that all parameters are the same).")
    private String energyRestart;

    /** -O or --noOriginal Do not include starting coordinates as their own rotamer. */
    @Option(names = {"-O",
        "--noOriginal"}, defaultValue = "false", description = "Do not include starting coordinates as their own rotamer.")
    private boolean noOriginal;

    /** -E or --decompose Print energy decomposition for the input structure (no optimization). */
    @Option(names = {"-E",
        "--decompose"}, defaultValue = "false", description = "Print energy decomposition for the input structure (no optimization!).")
    private boolean decompose;

    /**
     * --pH or --titrationPH Optimize the titration state of ASP, GLU, HIS and LYS residues.
     */
    @Option(names = {"--pH",
        "--titrationPH"}, paramLabel = "0", defaultValue = "0", description = " Optimize the titration state of ASP, GLU, HIS and LYS residues at the given pH (pH = 0 turns off titration")
    private double titrationPH;

    /**
     * --pH or --titrationPH Optimize the titration state of ASP, GLU, HIS and LYS residues.
     */
    @Option(names = {"--tR",
            "--titration"}, paramLabel = "false", defaultValue = "false", description = " Turn on titration state optimization")
    private boolean titration;

    /**
     * --mC or --monteCarlo Follow elimination criteria with 'n' Monte Carlo steps, or enumerate all
     * remaining conformations, whichever is smaller.
     */
    @Option(names = {"--mC",
        "--monteCarlo"}, paramLabel = "-1", defaultValue = "-1", description = "Follow elimination criteria with (n) Monte Carlo steps, or enumerate all remaining conformations, whichever is smaller.")
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

  }

  /**
   * Collection of ManyBody Box Optimization Options.
   */
  private static class BoxOptionGroup {

    /** -nB or --numBoxes Specify number of boxes along X, Y, and Z (default: '3,3,3'). */
    @Option(names = {"--nB",
        "--numBoxes"}, paramLabel = "3,3,3", defaultValue = "3,3,3", description = "Specify number of boxes along X, Y, and Z (default: 3,3,3)")
    private String numBoxes;

    /**
     * Result of parsing numBoxes flag.
     */
    private int[] numXYZBoxes;

    /** -bB or --boxBorderSize Extent of overlap between optimization boxes (default: 0.0 A). */
    @Option(names = {"--bB",
        "--boxBorderSize"}, paramLabel = "0.0", defaultValue = "0.0", description = "Extent of overlap between optimization boxes in Angstroms.")
    private double boxBorderSize;

    /**
     * -bL or --approxBoxLength Approximate side lengths of boxes to be constructed (over-rides
     * numXYZBoxes). Box sizes are rounded up to make a whole number of boxes along each axis
     * (default of 0 disables this function).
     */
    @Option(names = {"--bL",
        "--approxBoxLength"}, paramLabel = "20.0", defaultValue = "20.0", description = "Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes).")
    private double approxBoxLength;

    /**
     * -bC or --boxInclusionCriterion Criterion to use for adding residues to boxes. (1) uses C alpha
     * only (N1/9 for nucleic acids) (2) uses any atom. (3) uses any rotamer
     */
    @Option(names = {"--bC",
        "--boxInclusionCriterion"}, paramLabel = "1", defaultValue = "1", description = "Criterion to use for adding a residue to a box: (1) uses C alpha only (N1/9 for nucleic acids), (2) uses any atom, and (3) uses any rotamer")
    private int boxInclusionCriterion;

    /**
     * --iB or --initialBox Initial box to optimize.
     */
    @Option(names = {"--iB", "--initialBox"}, defaultValue = "0", description = "Initial box to optimize.")
    private int initialBox;

    /**
     * --fB or --boxFinal Final box to optimize.
     */
    @Option(names = {"--fB", "--finalBox"}, defaultValue = "2147483647", // Integer.MAX_VALUE
        description = "Final box to optimize.")
    private int finalBox;

    /**
     * --bT or --boxTitration Center boxes around titratable residues.
     */
    @Option(names = {"--bT", "--boxTitration"}, defaultValue = "false", // Integer.MAX_VALUE
            description = "Center boxes around titratable residues.")
    private boolean boxTitration;

  }

  /**
   * Collection of ManyBody Window Optimization Options.
   */
  private static class WindowOptionGroup {

    /** --window Size of the sliding window with respect to adjacent residues (default = 7). */
    @Option(names = {
        "--window"}, paramLabel = "7", defaultValue = "7", description = "Size of the sliding window with respect to adjacent residues.")
    private int window;

    /** --increment Sliding window increment (default = 3). */
    @Option(names = {
        "--increment"}, paramLabel = "3", defaultValue = "3", description = "Sliding window increment.")
    private int increment;

  }

  /**
   * Collection of ManyBody Energy Optimization Options.
   */
  private static class EnergyOptionGroup {

    /** --radius The sliding window and box cutoff radius (Angstroms). */
    @Option(names = {
        "--radius"}, paramLabel = "2.0", defaultValue = "2.0", description = "The sliding box and window cutoff radius (Angstroms).")
    private double cutoff;

    /** --tC or --twoBodyCutoff Cutoff distance for two-body interactions. */
    @Option(names = {"--tC",
        "--twoBodyCutoff"}, paramLabel = "3.0", defaultValue = "3.0", description = "Cutoff distance for two body interactions.")
    private double twoBodyCutoff;

    /** -T or --threeBody Include 3-Body interactions in the elimination criteria. */
    @Option(names = {"-T",
        "--threeBody"}, defaultValue = "false", description = "Include 3-Body interactions in the elimination criteria.")
    private boolean threeBody;

    /** --thC or --threeBodyCutoff Cutoff distance for three-body interactions. */
    @Option(names = {"--thC",
        "--threeBodyCutoff"}, paramLabel = "3.0", defaultValue = "3.0", description = "Cutoff distance for three-body interactions.")
    private double threeBodyCutoff;

    /** --pr or --prune Prune no clashes (0), only single clashes (1), or all clashes (2). */
    @Option(names = {"--pr",
        "--prune"}, paramLabel = "1", defaultValue = "1", description = "Prune no clashes (0), only single clashes (1), or all clashes (2)")
    private int prune;

    /**
     * --clashThreshold The threshold for pruning clashes. If two self-energies on the same residue
     * have an energy difference greater than 25 kcal/mol, the high energy rotamers get pruned.
     */
    @Option(names = {
        "--clashThreshold"}, paramLabel = "25.0", defaultValue = "25.0", description = "The threshold for pruning clashes.")
    private double clashThreshold;

    /**
     * --pairClashThreshold The threshold for pruning clashes. If two pair-energies on the same
     * residues have an energy difference greater than 25 kcal/mol, the high energy rotamers get
     * pruned.
     */
    @Option(names = {
        "--pairClashThreshold"}, paramLabel = "25.0", defaultValue = "25.0", description = "The threshold for pruning pair clashes.")
    private double pairClashThreshold;

    /** --radius The sliding window and box cutoff radius (Angstroms). */
    @Option(names = {
            "--kPH", "--pHRestraint"}, paramLabel = "0.0", defaultValue = "0.0", description = "Only allow titration state to change from" +
            "standard state is self energy exceeds the restraint.")
    private double pHRestraint = 0;
  }

  /**
   * Collection of ManyBody Residue Selection Options.
   */
  private static class ResidueOptionGroup {

    /** --ch or --chain Single character chain ID of the residues to optimize. */
    @Option(names = {"--ch",
        "--chain"}, paramLabel = "<A>", defaultValue = "-1", description = "Include only specified chain ID (default: all chains).")
    private String chain;

    /**
     * --sR or --start Starting residue to perform the optimization on.
     */
    @Option(names = {"--sR", "--start"}, defaultValue = "-2147483648", // Integer.MIN_VALUE
        description = "Starting residue to optimize (default: all residues).")
    private int start;

    /**
     * --fi or --final Final residue to perform the optimization.
     */
    @Option(names = {"--fR",
        "--final"}, paramLabel = "<final>", defaultValue = "2147483647", // Integer.MAX_VALUE
        description = "Final residue to optimize (default: all residues).")
    private int finish;

    /** --lR or --listResidues Choose a list of individual residues to optimize (eg. A11,A24,B40). */
    @Option(names = {"--lR",
        "--listResidues"}, paramLabel = "<list>", defaultValue = "none", description = "Select a list of residues to optimize (eg. A11,A24,B40).")
    private String listResidues;

    /** --oT or --onlyTitrtaion Rotamer optimize only titratable residues. */
    @Option(names = {"--oT",
            "--onlyTitration"}, paramLabel = "", defaultValue = "false", description = "Rotamer optimize only titratable residues.")
    private boolean onlyTitration;

    /** --iR or --interestedResidue Optimize rotamers within some distance of a specific residue. */
    @Option(names = {"--iR",
            "--interestedResidue"}, paramLabel = "", defaultValue = "-1", description = "Optimize rotamers within some distance of a specific residue.")
    private int interestedResidue = -1;

    /** --iC or --inclusionCutoff Distance which rotamers will be included when using only protons, titratable residues, or interested residue. */
    @Option(names = {"--iC",
            "--inclusionCutoff"}, paramLabel = "", defaultValue = "-1", description = "Distance which rotamers will be included when using only protons, titratable residues, or interested residue.")
    private double inclusionCutoff = -1;

  }

}
