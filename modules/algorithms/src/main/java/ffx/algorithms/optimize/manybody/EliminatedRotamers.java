// ******************************************************************************
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
// ******************************************************************************
package ffx.algorithms.optimize.manybody;

import ffx.algorithms.optimize.RotamerOptimization;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;

import static ffx.potential.bonded.Residue.ResidueType.NA;
import static java.lang.String.format;

public class EliminatedRotamers {

  private static final Logger logger = Logger.getLogger(EliminatedRotamers.class.getName());

  private final RotamerOptimization rO;
  private final DistanceMatrix dM;
  /**
   * A list of all residues being optimized. Note that Box and Window optimizations operate on
   * subsets of this list.
   */
  private final List<Residue> allResiduesList;
  /**
   * Maximum depth to check if a rotamer can be eliminated.
   */
  private final int maxRotCheckDepth;
  /**
   * Clash energy threshold (kcal/mole).
   */
  private final double clashThreshold;
  /**
   * Clash energy threshold (kcal/mole).
   */
  private final double pairClashThreshold;
  /**
   * Clash energy threshold (kcal/mol) for MultiResidues, which can have much more variation in self
   * and 2-Body energies.
   */
  private final double multiResClashThreshold;
  /**
   * Factor by which to multiply the pruning constraints for nucleic acids. nucleicPairsPruningFactor
   * is the arithmetic mean of 1.0 and the pruning factor, and is applied for AA-NA pairs.
   */
  private final double nucleicPruningFactor;

  private final double nucleicPairsPruningFactor;
  /**
   * Pair clash energy threshold (kcal/mol) for MultiResidues.
   */
  private final double multiResPairClashAddn;
  /**
   * Flag to prune clashes.
   */
  private final boolean pruneClashes;
  /**
   * Flag to prune pair clashes.
   */
  private final boolean prunePairClashes;
  /**
   * Flag to control the verbosity of printing.
   */
  private final boolean print;
  /**
   * Eliminated rotamers. [residue][rotamer]
   */
  public boolean[][] eliminatedSingles;
  /**
   * Eliminated rotamer pairs. [residue1][rotamer1][residue2][rotamer2]
   */
  public boolean[][][][] eliminatedPairs;
  /**
   * Pruned rotamers. Only for JUnit testing purposes.
   */
  public boolean[][] onlyPrunedSingles;
  /**
   * Pruned rotamer pairs. Only for JUnit testing purposes.
   */
  public boolean[][][][] onlyPrunedPairs;

  private EnergyExpansion eE;

  public EliminatedRotamers(RotamerOptimization rO, DistanceMatrix dM, List<Residue> allResiduesList,
                            int maxRotCheckDepth, double clashThreshold, double pairClashThreshold,
                            double multiResClashThreshold, double nucleicPruningFactor, double nucleicPairsPruningFactor,
                            double multiResPairClashAddn, boolean pruneClashes, boolean prunePairClashes, boolean print,
                            Residue[] residues) {
    this.rO = rO;
    this.dM = dM;
    this.allResiduesList = allResiduesList;
    this.maxRotCheckDepth = maxRotCheckDepth;
    this.clashThreshold = clashThreshold;
    this.pairClashThreshold = pairClashThreshold;
    this.multiResClashThreshold = multiResClashThreshold;
    this.nucleicPruningFactor = nucleicPruningFactor;
    this.nucleicPairsPruningFactor = nucleicPairsPruningFactor;
    this.multiResPairClashAddn = multiResPairClashAddn;
    this.pruneClashes = pruneClashes;
    this.prunePairClashes = prunePairClashes;
    this.print = print;
    allocateEliminationMemory(residues);
  }

  /**
   * Check for eliminated rotamer; true if eliminated.
   *
   * @param i  Residue i.
   * @param ri Rotamer ri.
   * @return True if rotamer eliminated.
   */
  public boolean check(int i, int ri) {
    if (eliminatedSingles == null) {
      return false;
    }
    return eliminatedSingles[i][ri];
  }

  /**
   * Check for eliminated rotamer pair; true if eliminated.
   *
   * @param i  Residue i.
   * @param ri Rotamer ri.
   * @param j  Residue j.
   * @param rj Rotamer rj.
   * @return True if eliminated pair.
   */
  public boolean check(int i, int ri, int j, int rj) {
    if (eliminatedPairs == null) {
      return false;
    }
    // If j is an earlier residue than i, swap j with i, as eliminated
    // rotamers are stored with the earlier residue listed first.
    if (j < i) {
      int ii = i;
      int iri = ri;
      i = j;
      ri = rj;
      j = ii;
      rj = iri;
    }
    return eliminatedPairs[i][ri][j][rj];
  }

  /**
   * Check for pruned rotamer pair; true if eliminated. Only used during testing.
   *
   * @param i  Residue i.
   * @param ri Rotamer ri.
   * @param j  Residue j.
   * @param rj Rotamer rj.
   * @return a boolean.
   */
  public boolean checkPrunedPairs(int i, int ri, int j, int rj) {
    if (onlyPrunedPairs == null) {
      return false;
    }

    if (j < i) {
      int ii = i;
      int iri = ri;
      i = j;
      ri = rj;
      j = ii;
      rj = iri;
    }
    return onlyPrunedPairs[i][ri][j][rj];
  }

  /**
   * Check for pruned rotamer; true if eliminated. Only used during testing.
   *
   * @param i  The residue.
   * @param ri The rotamer.
   * @return a boolean.
   */
  public boolean checkPrunedSingles(int i, int ri) {
    if (onlyPrunedSingles == null) {
      return false;
    }
    return onlyPrunedSingles[i][ri];
  }

  /**
   * Checks to see if any eliminations with j,rj have occurred; assumes i,ri self has already been
   * checked. Checks j,rj self and i,ri,j,rj 2-Body. The intent is to be part of a loop over
   * i,ri,j,rj, and check for eliminations at the j,rj point.
   *
   * @param i  Residue i
   * @param ri Rotamer ri
   * @param j  Residue j
   * @param rj Rotamer rj
   * @return j eliminated with i
   */
  public boolean checkToJ(int i, int ri, int j, int rj) {
    return (check(j, rj) || check(i, ri, j, rj));
  }

  /**
   * Checks to see if any eliminations with k,rk have occurred; assumes i,ri,j,rj 2-Body has already
   * been checked. Checks the k,rk self, all pairs with k,rk, and the i,ri,j,rj,k,rk 3-Body. The
   * intent is to be part of a loop over i,ri,j,rj,k,rk, and check for eliminations at the k,rk
   * point.
   *
   * @param i  Residue i
   * @param ri Rotamer ri
   * @param j  Residue j
   * @param rj Rotamer rj
   * @param k  Residue k
   * @param rk Rotamer rk
   * @return k eliminated with i,j
   */
  public boolean checkToK(int i, int ri, int j, int rj, int k, int rk) {
    return (check(k, rk) || check(i, ri, k, rk) || check(j, rj, k, rk));
  }

  /**
   * Checks to see if any eliminations with l,rl have occurred; assumes i,ri,j,rj,k,rk 3-Body has
   * already been checked. Checks the l,rl self, all pairs with l,rl, all triples with l,rl, and the
   * 4-Body. The intent is to be part of a loop over i,ri,j,rj,k,rk,l,rl, and check for eliminations
   * at the l,rl point.
   *
   * @param i  Residue i
   * @param ri Rotamer ri
   * @param j  Residue j
   * @param rj Rotamer rj
   * @param k  Residue k
   * @param rk Rotamer rk
   * @param l  Residue l
   * @param rl Rotamer rl
   * @return l eliminated with i,j,k
   */
  public boolean checkToL(int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
    return (check(l, rl) || check(i, ri, l, rl) || check(j, rj, l, rl) || check(k, rk, l, rl));
  }

  /**
   * Safe method to eliminate a rotamer: will not eliminate if there are no alternate rotamers for
   * residue i, or if i-ri is already eliminated.
   *
   * @param residues Residues under consideration.
   * @param i        A residue index based on the current residue list.
   * @param ri       A rotamer to attempt elimination of.
   * @param verbose  Request verbose logging.
   * @return If the rotamer was eliminated.
   */
  public boolean eliminateRotamer(Residue[] residues, int i, int ri, boolean verbose) {
    // Check if rotamer (i, ri) has already been eliminated.
    if (check(i, ri)) {
      return false;
    }

    // Make sure at least one rotamer rii != ri is left.
    int[] validRots = rotamerCount(residues, i);
    int rotCount = 0;
    for (int rii = 0; rii < validRots.length; rii++) {
      if (rii != ri) {
        ++rotCount;
      }
    }

    if (rotCount == 0) {
      // No valid rotamers other than ri are left!
      return false;
    }

    eliminatedSingles[i][ri] = true;

    if (verbose) {
      Rotamer[] rotamers = residues[i].getRotamers();
      rO.logIfRank0(
          format(" Rotamer (%8s,%2d) eliminated (%2d left).", residues[i].toString(rotamers[ri]), ri,
              rotCount));
    }
    int eliminatedPairs = eliminateRotamerPairs(residues, i, ri, verbose);
    if (eliminatedPairs > 0 && verbose) {
      rO.logIfRank0(format("  Eliminated %2d rotamer pairs.", eliminatedPairs));
    }
    return true;
  }

  public boolean eliminateRotamerPair(Residue[] residues, int i, int ri, int j, int rj,
                                      boolean verbose) {
    if (i > j) {
      int ii = i;
      int iri = ri;
      i = j;
      ri = rj;
      j = ii;
      rj = iri;
    }

    if (!check(i, ri, j, rj)) {
      eliminatedPairs[i][ri][j][rj] = true;
      if (verbose) {
        Rotamer[] rotI = residues[i].getRotamers();
        Rotamer[] rotJ = residues[j].getRotamers();
        rO.logIfRank0(format("  Rotamer pair eliminated: [(%8s,%2d) (%8s,%2d)]",
            residues[i].toString(rotI[ri]), ri, residues[j].toString(rotJ[rj]), rj));
      }
      return true;
    } else {
      return false;
    }
  }

  public int eliminateRotamerPairs(Residue[] residues, int i, int ri, boolean verbose) {
    int eliminatedPairs = 0;
    for (int j = 0; j < residues.length; j++) {
      if (j == i) {
        continue;
      }
      Residue resJ = residues[j];
      int lenRj = resJ.getRotamers().length;
      for (int rj = 0; rj < lenRj; rj++) {
        if (eliminateRotamerPair(residues, i, ri, j, rj, verbose)) {
          ++eliminatedPairs;
        }
      }
    }
    return eliminatedPairs;
  }

  /**
   * Method to check if pairs elimination for some residue pair has enabled a singles rotamer
   * elimination by eliminating all ri-rj for some ri or some rj.
   *
   * @param residues Residues under consideration.
   * @param i        A residue index.
   * @param j        A residue index j!=i
   * @return If any singletons were eliminated.
   */
  public boolean pairsToSingleElimination(Residue[] residues, int i, int j) {
    assert i != j;
    assert i < residues.length;
    assert j < residues.length;

    Residue residueI = residues[i];
    Residue residueJ = residues[j];
    Rotamer[] rotI = residueI.getRotamers();
    Rotamer[] rotJ = residueJ.getRotamers();
    int lenRi = rotI.length;
    int lenRj = rotJ.length;
    boolean eliminated = false;

    // Now check ris with no remaining pairs to j.
    for (int ri = 0; ri < lenRi; ri++) {
      if (check(i, ri)) {
        continue;
      }
      boolean pairRemaining = false;
      for (int rj = 0; rj < lenRj; rj++) {
        if (!check(j, rj) && !check(i, ri, j, rj)) {
          pairRemaining = true;
          break;
        }
      }
      if (!pairRemaining) {
        if (eliminateRotamer(residues, i, ri, print)) {
          eliminated = true;
          rO.logIfRank0(format(" Eliminating rotamer %s-%d with no remaining pairs to residue %s.",
              residueI.toString(rotI[ri]), ri, residueJ));
        } else {
          rO.logIfRank0(
              format(" Already eliminated rotamer %s-%d with no remaining pairs to residue %s.",
                  residueI.toString(rotI[ri]), ri, residueJ), Level.WARNING);
        }
      }
    }

    // Check residue j rotamers with no remaining pairs to residue i.
    for (int rj = 0; rj < lenRj; rj++) {
      if (check(j, rj)) {
        continue;
      }
      boolean pairRemaining = false;
      for (int ri = 0; ri < lenRi; ri++) {
        if (!check(i, ri) && !check(i, ri, j, rj)) {
          pairRemaining = true;
          break;
        }
      }
      if (!pairRemaining) {
        if (eliminateRotamer(residues, j, rj, print)) {
          eliminated = true;
          rO.logIfRank0(format(" Eliminating rotamer %s-%d with no remaining pairs to residue %s.",
              residueJ.toString(rotJ[rj]), rj, residueI));
        } else {
          rO.logIfRank0(
              format(" Already eliminated rotamer J %s-%d with no remaining pairs to residue %s.",
                  residueJ.toString(rotJ[rj]), rj, residueI), Level.WARNING);
        }
      }
    }

    return eliminated;
  }

  /**
   * Pre-prunes any pairs that have a pair-energy of Double.NaN before pruning and eliminations
   * happen.
   *
   * @param residues Array of all residues.
   */
  public void prePrunePairs(Residue[] residues) {
    // Pre-Prune if pair-energy is Double.NaN.
    // Loop over first residue.
    int nResidues = residues.length;
    for (int i = 0; i < nResidues - 1; i++) {
      Residue resi = residues[i];
      Rotamer[] rotI = resi.getRotamers();
      int ni = rotI.length;
      // Loop over second residue.
      for (int j = i + 1; j < nResidues; j++) {
        Residue resJ = residues[j];
        Rotamer[] rotJ = resJ.getRotamers();
        int nj = rotJ.length;
        // Loop over the rotamers for residue i.
        for (int ri = 0; ri < ni; ri++) {
          if (!validRotamer(residues, i, ri)) {
            continue;
          }
          // Loop over rotamers for residue j.
          for (int rj = 0; rj < nj; rj++) {
            if (!validRotamer(residues, j, rj) || check(i, ri, j, rj)) {
              continue;
            }
            if (!check(i, ri, j, rj) && Double.isNaN(eE.get2Body(i, ri, j, rj))) {
              rO.logIfRank0(format(
                  " Rotamer Pair (%7s,%2d) (%7s,%2d) 2-body energy %12.4f pre-pruned since energy is NaN.",
                  i, ri, j, rj, eE.get2Body(i, ri, j, rj)));
              eliminateRotamerPair(residues, i, ri, j, rj, print);
            }
          }
        }
      }
    }
  }

  /**
   * Pre-prunes any selves that have a self-energy of Double.NaN before pruning and eliminations
   * happen.
   *
   * @param residues Array of all residues.
   */
  public void prePruneSelves(Residue[] residues) {
    // Pre-Prune if self-energy is Double.NaN.
    for (int i = 0; i < residues.length; i++) {
      Residue residue = residues[i];
      Rotamer[] rotamers = residue.getRotamers();
      int nRot = rotamers.length;
      for (int ri = 0; ri < nRot; ri++) {
        if (!check(i, ri) && Double.isNaN(eE.getSelf(i, ri))) {
          rO.logIfRank0(
              format(" Rotamer (%7s,%2d) self-energy %12.4f pre-pruned since energy is NaN.",
                  residue, ri, eE.getSelf(i, ri)));
          eliminateRotamer(residues, i, ri, false);
        }
      }
    }
  }

  /**
   * Prunes rotamer ri of residue i if all ri-j pair energies are worse than the best i-j pair by
   * some threshold value; additionally prunes ri-rj pairs if they exceed the best i-j pair by a
   * greater threshold value; additionally performs this in reverse (searches over j-i).
   *
   * @param residues Residues whose rotamers are to be pruned.
   */
  public void prunePairClashes(Residue[] residues) {
    if (!prunePairClashes) {
      return;
    }
    int nResidues = residues.length;

    for (int i = 0; i < nResidues - 1; i++) {
      Residue residueI = residues[i];
      Rotamer[] rotI = residueI.getRotamers();
      int lenRi = rotI.length;
      int indI = allResiduesList.indexOf(residueI);
      for (int j = i + 1; j < nResidues; j++) {
        Residue residueJ = residues[j];
        Rotamer[] rotJ = residueJ.getRotamers();
        int lenRj = rotJ.length;
        int indJ = allResiduesList.indexOf(residueJ);

        double minPair = Double.MAX_VALUE;
        int minRI = -1;
        int minRJ = -1;

        boolean cutoffPair = true;
        for (int ri = 0; ri < lenRi; ri++) {
          if (check(i, ri)) {
            continue;
          }
          for (int rj = 0; rj < lenRj; rj++) {
            if (check(j, rj) || check(i, ri, j, rj) || dM.checkPairDistThreshold(indI, ri, indJ,
                rj)) {
              continue;
            }
            cutoffPair = false;
            double pairEnergy = eE.get2Body(i, ri, j, rj) + eE.getSelf(i, ri) + eE.getSelf(j, rj);
            assert Double.isFinite(pairEnergy);
            if (pairEnergy < minPair) {
              minPair = pairEnergy;
              minRI = ri;
              minRJ = rj;
            }
          }
        }
        if (cutoffPair) {
          // Under this branch: no rotamer pairs that clear the distance threshold.
          continue;
        }
        assert (minRI >= 0 && minRJ >= 0); // Otherwise, i and j are not on a well-packed backbone.

        // Calculate all the modifiers to the pair clash elimination threshold.
        double threshold = pairClashThreshold;
        if (residueI instanceof MultiResidue) {
          threshold += multiResPairClashAddn;
        }
        if (residueJ instanceof MultiResidue) {
          threshold += multiResPairClashAddn;
        }
        int numNARes =
            (residueI.getResidueType() == NA ? 1 : 0) + (residueJ.getResidueType() == NA ? 1 : 0);
        switch (numNARes) {
          case 0:
            break;
          case 1:
            threshold *= nucleicPairsPruningFactor;
            break;
          case 2:
            threshold *= nucleicPruningFactor;
            break;
          default:
            throw new ArithmeticException(" RotamerOptimization.prunePairClashes() has somehow "
                + "found less than zero or more than two nucleic acid residues in a pair of"
                + " residues. This result should be impossible.");
        }
        double toEliminate = threshold + minPair;

        for (int ri = 0; ri < lenRi; ri++) {
          if (check(i, ri)) {
            continue;
          }
          for (int rj = 0; rj < lenRj; rj++) {
            if (check(j, rj) || check(i, ri, j, rj)) {
              continue;
            }
            double pairEnergy = eE.get2Body(i, ri, j, rj) + eE.getSelf(i, ri) + eE.getSelf(j, rj);
            assert Double.isFinite(pairEnergy);
            if (pairEnergy > toEliminate) {
              rO.logIfRank0(
                  format(" Pruning pair %s-%d %s-%d by %s-%d %s-%d; energy %s > " + "%s + %s",
                      residueI.toString(rotI[ri]), ri, residueJ.toString(rotJ[rj]), rj,
                      residueI.toString(rotI[minRI]), minRI, residueJ.toString(rotJ[minRJ]), minRJ,
                      rO.formatEnergy(pairEnergy), rO.formatEnergy(threshold),
                      rO.formatEnergy(minPair)));
            }
          }
        }

        pairsToSingleElimination(residues, i, j);
      }
    }
  }

  /**
   * Uses calculated energies to prune rotamers based on a threshold distance from that residue's
   * minimum energy rotamer (by default 20 kcal/mol). The threshold can be modulated by presence of
   * nucleic acids or MultiResidues, which require more generous pruning criteria.
   *
   * @param residues Residues to prune rotamers over.
   */
  public void pruneSingleClashes(Residue[] residues) {
    if (!pruneClashes) {
      return;
    }
    for (int i = 0; i < residues.length; i++) {
      Residue residue = residues[i];
      Rotamer[] rotamers = residue.getRotamers();
      int nRot = rotamers.length;
      double minEnergy = Double.MAX_VALUE;
      int minRot = -1;
      for (int ri = 0; ri < nRot; ri++) {
        if (!check(i, ri) && eE.getSelf(i, ri) < minEnergy) {
          minEnergy = eE.getSelf(i, ri);
          minRot = ri;
        }
      }

      double energyToPrune =
          (residue instanceof MultiResidue) ? multiResClashThreshold : clashThreshold;
      energyToPrune =
          (residue.getResidueType() == NA) ? energyToPrune * nucleicPruningFactor : energyToPrune;
      energyToPrune += minEnergy;

      for (int ri = 0; ri < nRot; ri++) {
        if (!check(i, ri) && (eE.getSelf(i, ri) > energyToPrune)) {
          if (eliminateRotamer(residues, i, ri, print)) {
            rO.logIfRank0(format("  Rotamer (%7s,%2d) self-energy %s pruned by (%7s,%2d) %s.",
                residue.toString(rotamers[ri]), ri, rO.formatEnergy(eE.getSelf(i, ri)),
                residue.toString(rotamers[minRot]), minRot, rO.formatEnergy(minEnergy)));
          }
        }
      }
    }
  }

  public void setEnergyExpansion(EnergyExpansion eE) {
    this.eE = eE;
  }

  @Override
  public String toString() {
    int rotamerCount = 0;
    int pairCount = 0;
    int singles = 0;
    int pairs = 0;
    int nRes = eliminatedSingles.length;
    for (int i = 0; i < nRes; i++) {
      int nRotI = eliminatedSingles[i].length;
      rotamerCount += nRotI;
      for (int ri = 0; ri < nRotI; ri++) {
        if (eliminatedSingles[i][ri]) {
          singles++;
        }
        for (int j = i + 1; j < nRes; j++) {
          int nRotJ = eliminatedPairs[i][ri][j].length;
          pairCount += nRotJ;
          for (int rj = 0; rj < nRotJ; rj++) {
            if (eliminatedPairs[i][ri][j][rj]) {
              pairs++;
            }
          }
        }
      }
    }
    return format(" %d out of %d rotamers eliminated.\n", singles, rotamerCount) + format(
        " %d out of %d rotamer pairs eliminated.", pairs, pairCount);
  }

  public boolean validateDEE(Residue[] residues) {
    int nRes = eliminatedSingles.length;
    // Validate residues
    for (int i = 0; i < nRes; i++) {
      Residue residueI = residues[i];
      int ni = eliminatedSingles[i].length;
      boolean valid = false;
      for (int ri = 0; ri < ni; ri++) {
        if (!check(i, ri)) {
          valid = true;
        }
      }
      if (!valid) {
        logger.severe(
            format(" Coding error: all %d rotamers for residue %s eliminated.", ni, residueI));
      }
    }

    // Validate pairs
    for (int i = 0; i < nRes; i++) {
      Residue residueI = residues[i];
      Rotamer[] rotI = residueI.getRotamers();
      int ni = rotI.length;
      for (int j = i + 1; j < nRes; j++) {
        Residue residueJ = residues[j];
        Rotamer[] rotJ = residueJ.getRotamers();
        int nj = rotJ.length;
        boolean valid = false;
        for (int ri = 0; ri < ni; ri++) {
          for (int rj = 0; rj < nj; rj++) {
            if (!check(i, ri, j, rj)) {
              valid = true;
            }
          }
        }
        if (!valid) {
          logger.severe(format(" Coding error: all pairs for %s with residue %s eliminated.",
              residueI.toFormattedString(false, true), residueJ));
        }
      }
    }

    return true;
  }

  /**
   * allocateEliminationMemory.
   *
   * @param residues an array of {@link ffx.potential.bonded.Residue} objects.
   */
  private void allocateEliminationMemory(Residue[] residues) {
    int nRes = residues.length;
    eliminatedSingles = new boolean[nRes][];
    eliminatedPairs = new boolean[nRes][][][];
    // Loop over residues.
    rO.logIfRank0("\n     Residue  Nrot");
    for (int i = 0; i < nRes; i++) {
      Residue residueI = residues[i];
      Rotamer[] rotamersI = residueI.getRotamers();
      int lenRi = rotamersI.length;
      rO.logIfRank0(format(" %3d %8s %4d", i + 1, residueI.toFormattedString(false, true), lenRi));
      eliminatedSingles[i] = new boolean[lenRi];
      eliminatedPairs[i] = new boolean[lenRi][][];
      // Loop over the set of rotamers for residue i.
      for (int ri = 0; ri < lenRi; ri++) {
        eliminatedSingles[i][ri] = false;
        eliminatedPairs[i][ri] = new boolean[nRes][];
        for (int j = i + 1; j < nRes; j++) {
          Residue residueJ = residues[j];
          Rotamer[] rotamersJ = residueJ.getRotamers();
          int lenRj = rotamersJ.length;
          eliminatedPairs[i][ri][j] = new boolean[lenRj];
          for (int rj = 0; rj < lenRj; rj++) {
            eliminatedPairs[i][ri][j][rj] = false;
          }
        }
      }
    }
  }

  /**
   * Validate residue i with rotamer ri.
   *
   * @param residues The residues being optimized.
   * @param i        The residue to validate.
   * @param ri       The rotamer to validate.
   * @return The status of this rotamer.
   */
  private boolean validRotamer(Residue[] residues, int i, int ri) {
    // Return false if this rotamer has been eliminated.
    if (check(i, ri)) {
      return false;
    }

    if (maxRotCheckDepth > 1) {
      // Loop over all residues to check for valid pairs and triples.
      int n = residues.length;
      for (int j = 0; j < n; j++) {
        if (j == i) {
          continue;
        }
        if (rotamerPairCount(residues, i, ri, j) == 0) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Count the rotamers remaining for residue i.
   *
   * @param residues Residue array.
   * @param i        The residue number to examine.
   * @return The remaining valid rotamers.
   */
  private int[] rotamerCount(Residue[] residues, int i) {
    int nRes = residues.length;
    Rotamer[] rotI = residues[i].getRotamers();
    int ni = rotI.length;

    if (maxRotCheckDepth == 0) {
      // Short-circuit on all its rotamers.
      return IntStream.range(0, ni).toArray();
    }

    return IntStream.range(0, ni).filter((int ri) -> {
      if (check(i, ri)) {
        return false;
      }
      if (maxRotCheckDepth > 1) {
        // Check that rotamer ri has valid pairs with all other residues.
        for (int j = 0; j < nRes; j++) {
          if (i == j) {
            continue;
          }
          if (rotamerPairCount(residues, i, ri, j) == 0) {
            return false;
          }
        }
      }
      return true;
    }).toArray();
  }

  /**
   * Validate rotamer pair (i, ri) and (j, rj).
   *
   * @param residues The residues being optimized.
   * @param i        The first residue to validate.
   * @param ri       The first rotamer to validate.
   * @param j        The 2nd residue to validate.
   * @param rj       The 2nd rotamer to validate.
   * @return The status of this rotamer.
   */
  private boolean validRotamerPair(Residue[] residues, int i, int ri, int j, int rj) {
    // Residues i and j must be different.
    if (i == j) {
      return false;
    }

    // Return false if either rotamer is not valid.
    if (!validRotamer(residues, i, ri) || !validRotamer(residues, j, rj)) {
      return false;
    }

    // Return false if the rotamer pair has been eliminated.
    if (check(i, ri, j, rj)) {
      return false;
    }

    if (maxRotCheckDepth > 1) {
      // Loop over all residues to check for valid triples.
      int n = residues.length;
      for (int k = 0; k < n; k++) {
        if (k == i || k == j) {
          continue;
        }
        // There must be at least one valid rotamer triple between (i, ri) and (j, rj) with residue
        // k.
        if (rotamerTripleCount(residues, i, ri, j, rj, k) == 0) {
          // Eliminate the pair?
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Count the rotamer pairs remaining for (residue i, rotamer ri) and residue j.
   *
   * @param residues Residue array.
   * @param i        The first residue to examine.
   * @param ri       The rotamer for the first residue.
   * @param j        The second residue to examine.
   * @return The remaining rotamer pair count.
   */
  private int rotamerPairCount(Residue[] residues, int i, int ri, int j) {
    if (i == j || check(i, ri)) {
      return 0;
    }
    int pairCount = 0;
    Rotamer[] rotJ = residues[j].getRotamers();
    int nj = rotJ.length;
    // Loop over all rotamers for residue j.
    for (int rj = 0; rj < nj; rj++) {
      // Check for a valid rotamer pair.
      if (!check(j, rj) && !check(i, ri, j, rj)) {
        // Loop over all residues k to check for valid rotamer triples.
        int nRes = residues.length;
        boolean valid = true;
        if (maxRotCheckDepth > 2) {
          for (int k = 0; k < nRes; k++) {
            if (k == i || k == j) {
              continue;
            }
            if (rotamerTripleCount(residues, i, ri, j, rj, k) == 0) {
              valid = false;
            }
          }
        }
        if (valid) {
          pairCount++;
        }
      }
    }
    return pairCount;
  }

  /**
   * Count the rotamer triples remaining for (residue i, rotamer ri) and (residue j, rotamer rj) with
   * residue k.
   *
   * @param residues Residue array.
   * @param i        The first residue to examine.
   * @param ri       The rotamer for the first residue.
   * @param j        The second residue to examine.
   * @param rj       The rotamer for the first residue.
   * @param k        The third residue.
   * @return The remaining rotamer triples count.
   */
  private int rotamerTripleCount(Residue[] residues, int i, int ri, int j, int rj, int k) {
    if (i == j || i == k || j == k) {
      return 0;
    }
    int tripleCount = 0;
    Rotamer[] rotK = residues[k].getRotamers();
    int nk = rotK.length;
    // Check that each rotamer and their pair have not been eliminated.
    if (!check(i, ri) && !check(j, rj) && !check(i, ri, j, rj)) {
      // Loop over all rotamers for residue k.
      for (int rk = 0; rk < nk; rk++) {
        // Check for a valid rotamer triple.
        if (!check(k, rk) && !check(i, ri, k, rk) && !check(j, rj, k, rk)) {
          tripleCount++;
        }
      }
    }
    return tripleCount;
  }
}
