// ******************************************************************************
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
// ******************************************************************************
package ffx.algorithms.optimize.manybody;

import static java.lang.String.format;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.ParallelTeam;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.optimize.RotamerOptimization;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.nonbonded.NeighborList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

/**
 * Calculates a residue-residue distance matrix.
 *
 * <p>Residue-residue distance is defined as the shortest atom-atom distance in any possible
 * rotamer-rotamer pair if the residues are neighbors (central atom-central atom distances are
 * within a cutoff). Otherewise, distances are set to a default of Double.MAX_VALUE.
 *
 * <p>The intent of using a neighbor list is to avoid tediously searching rotamer- rotamer pairs
 * when two residues are so far apart we will never need the exact distance. We use the distance
 * matrix for adding residues to the sliding window and determining whether to set 3-body energy to
 * 0.0.
 *
 * <p>If the central atoms are too distant from each other, we can safely assume no atoms will ever
 * be close enough for addition to sliding window or to cause explicit calculation of 3-body energy.
 */
public class DistanceMatrix {

  private static final Logger logger = Logger.getLogger(DistanceMatrix.class.getName());

  private final RotamerOptimization rO;
  /** MolecularAssembly to perform rotamer optimization on. */
  private final MolecularAssembly molecularAssembly;
  /** AlgorithmListener who should receive updates as the optimization runs. */
  private final AlgorithmListener algorithmListener;
  /**
   * An array of all residues being optimized. Note that Box and Window optimizations operate on
   * subsets of this list.
   */
  private final Residue[] allResiduesArray;
  /**
   * A list of all residues being optimized. Note that Box and Window optimizations operate on
   * subsets of this list.
   */
  private final List<Residue> allResiduesList;
  /** Number of residues being optimized. */
  private final int numResidues;

  private final RotamerLibrary library;
  /** Default distance method is to find the shortest distance between residues. */
  private final RotamerOptimization.DistanceMethod distanceMethod;
  /** The distance that the distance matrix checks for. */
  private final double distance;
  /** Two-Body cutoff distance. */
  private final double twoBodyCutoffDist;
  /** Three-body cutoff distance. */
  private final double threeBodyCutoffDist;
  /**
   * Flag to load the distance matrix as needed; if false, matrix is prefilled at the beginning of
   * rotamer optimization.
   */
  private final boolean lazyMatrix;
  /** Flag to indicate use of forced residues during the sliding window algorithm. */
  private final boolean useForcedResidues;
  /**
   * The minimum distance between atoms of a residue pair, taking into account interactions with
   * symmetry mates.
   *
   * <p>[residue1][rotamer1][residue2][rotamer2]
   */
  private double[][][][] distanceMatrix;

  public DistanceMatrix(
      RotamerOptimization rO,
      MolecularAssembly molecularAssembly,
      AlgorithmListener algorithmListener,
      Residue[] allResiduesArray,
      List<Residue> allResiduesList,
      RotamerLibrary library,
      RotamerOptimization.DistanceMethod distanceMethod,
      double distance,
      double twoBodyCutoffDist,
      double threeBodyCutoffDist,
      boolean lazyMatrix,
      boolean useForcedResidues) {
    this.rO = rO;
    this.molecularAssembly = molecularAssembly;
    this.algorithmListener = algorithmListener;
    this.allResiduesArray = allResiduesArray;
    this.allResiduesList = allResiduesList;
    this.numResidues = allResiduesArray.length;
    this.library = library;
    this.distanceMethod = distanceMethod;
    this.distance = distance;
    this.twoBodyCutoffDist = twoBodyCutoffDist;
    this.threeBodyCutoffDist = threeBodyCutoffDist;
    this.lazyMatrix = lazyMatrix;
    this.useForcedResidues = useForcedResidues;

    distanceMatrix();
  }

  /**
   * Returns true if all values (usually distances) are both finite and not set to Double.MAX_VALUE.
   *
   * @param values Values to check.
   * @return If all values are regular, finite values.
   */
  private static boolean areFiniteAndNotMax(double... values) {
    return Arrays.stream(values)
        .allMatch((double val) -> Double.isFinite(val) && val < Double.MAX_VALUE);
  }

  /**
   * Gets the raw distance between two rotamers using lazy loading of the distance matrix. Intended
   * uses: helper method for get2BodyDistance, and for checking for superpositions.
   *
   * @param i A residue index.
   * @param ri A rotamer index for i.
   * @param j A residue index j!=i.
   * @param rj A rotamer index for j.
   * @return dist(i, ri, j, rj), ignoring distance matrix method used.
   */
  public double checkDistMatrix(int i, int ri, int j, int rj) {
    if (i > j) {
      int temp = i;
      i = j;
      j = temp;
      temp = ri;
      ri = rj;
      rj = temp;
    }
    double dist = distanceMatrix[i][ri][j][rj];
    if (dist < 0) {
      dist = evaluateDistance(i, ri, j, rj);
      distanceMatrix[i][ri][j][rj] = dist;
    }
    return dist;
  }

  /**
   * Checks if the i,ri,j,rj pair exceeds the pair distance thresholds.
   *
   * @param i  A residue index.
   * @param ri A rotamer index for i.
   * @param j  A residue index j!=i.
   * @param rj A rotamer index for j.
   * @return If i,ri,j,rj is greater than the threshold distance.
   */
  public boolean checkPairDistThreshold(int i, int ri, int j, int rj) {
    if (twoBodyCutoffDist <= 0 || !Double.isFinite(twoBodyCutoffDist)) {
      return false;
    }
    return (get2BodyDistance(i, ri, j, rj) > twoBodyCutoffDist);
  }

  /**
   * Checks if the i,ri,j,rj,k,rk,l,rl quad exceeds the 3-body threshold, or if any component
   * exceeds the pair/triple distance thresholds.
   *
   * @param i A residue index.
   * @param ri A rotamer index for i.
   * @param j A residue index j!=i.
   * @param rj A rotamer index for j.
   * @param k A residue index k!=i, k!=j.
   * @param rk A rotamer index for k.
   * @param l A residue index l!=i, l!=j, l!=k.
   * @param rl A rotamer index for l.
   * @return If i,ri,j,rj,k,rk,l,rl is greater than the threshold distances.
   */
  public boolean checkQuadDistThreshold(
      int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
    if (checkTriDistThreshold(i, ri, j, rj, k, rk)
        || checkTriDistThreshold(i, ri, j, rj, l, rl)
        || checkTriDistThreshold(i, ri, k, rk, l, rl)
        || checkTriDistThreshold(j, rj, k, rk, l, rl)) {
      return true;
    }
    // Use the 3-body cutoff distance for now.
    if (threeBodyCutoffDist <= 0 || !Double.isFinite(threeBodyCutoffDist)) {
      return false;
    }
    return get4BodyDistance(i, ri, j, rj, k, rk, l, rl) > threeBodyCutoffDist;
  }

  /**
   * Checks if the i,ri,j,rj,k,rk triple exceeds the 3-body threshold, or if any component exceeds
   * the pair distance threshold.
   *
   * @param i A residue index.
   * @param ri A rotamer index for i.
   * @param j A residue index j!=i.
   * @param rj A rotamer index for j.
   * @param k A residue index k!=i, k!=j.
   * @param rk A rotamer index for k.
   * @return If i,ri,j,rj,k,rk is greater than the threshold distances.
   */
  public boolean checkTriDistThreshold(int i, int ri, int j, int rj, int k, int rk) {
    if (checkPairDistThreshold(i, ri, j, rj)
        || checkPairDistThreshold(i, ri, k, rk)
        || checkPairDistThreshold(j, rj, k, rk)) {
      return true;
    }
    if (threeBodyCutoffDist <= 0 || !Double.isFinite(threeBodyCutoffDist)) {
      return false;
    }
    return (get3BodyDistance(i, ri, j, rj, k, rk) > threeBodyCutoffDist);
  }

  /**
   * Checks the distance matrix, finding the shortest distance between two residues' rotamers or any
   * rotamers from two residues under any symmetry operator; will evaluate this if distance matrix
   * not already filled. Default is to find the shortest distance between two residues rather than
   * between two rotamers.
   *
   * @param i Residue i
   * @param ri Rotamer for i
   * @param j Residue j
   * @param rj Rotamer for j
   * @return Shortest distance
   */
  public double get2BodyDistance(int i, int ri, int j, int rj) {
    if (i > j) {
      int temp = i;
      i = j;
      j = temp;
      temp = ri;
      ri = rj;
      rj = temp;
    }

    if (distanceMethod == RotamerOptimization.DistanceMethod.ROTAMER) {
      return checkDistMatrix(i, ri, j, rj);
    } else {
      return getResidueDistance(i, ri, j, rj);
    }
  }

  /**
   * Returns the RMS separation distance for the closest rotamers of three residues' 2-body
   * distances. Defaults to Double.MAX_VALUE when there are pair distances outside cutoffs.
   *
   * @param i Residue i
   * @param ri Rotamer for i
   * @param j Residue j
   * @param rj Rotamer for j
   * @param k Residue k
   * @param rk Rotamer for k
   * @return RMS separation distance
   */
  public double get3BodyResidueDistance(int i, int ri, int j, int rj, int k, int rk) {
    double ij = getResidueDistance(i, ri, j, rj);
    double ik = getResidueDistance(i, ri, k, rk);
    double jk = getResidueDistance(j, rj, k, rk);
    if (areFiniteAndNotMax(ij, ik, jk)) {
      return sqrt((ij * ij + ik * ik + jk * jk) / 3.0);
    } else {
      return Double.MAX_VALUE;
    }
  }

  /**
   * Returns the RMS separation distance for the closest rotamers of 6 2-body distances from four
   * residues. Defaults to Double.MAX_VALUE when there are pair distances outside cutoffs.
   *
   * @param i Residue i
   * @param ri Rotamer for i
   * @param j Residue j
   * @param rj Rotamer for j
   * @param k Residue k
   * @param rk Rotamer for k
   * @param l Residue l
   * @param rl Rotamer for l
   * @return RMS separation distance
   */
  public double get4BodyResidueDistance(
      int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
    double ij = getResidueDistance(i, ri, j, rj);
    double ik = getResidueDistance(i, ri, k, rk);
    double il = getResidueDistance(i, ri, l, rl);
    double jk = getResidueDistance(j, rj, k, rk);
    double jl = getResidueDistance(j, rj, l, rl);
    double kl = getResidueDistance(k, rk, l, rl);
    if (areFiniteAndNotMax(ij, ik, il, jk, jl, kl)) {
      return sqrt((ij * ij + ik * ik + il * il + jk * jk + jl * jl + kl * kl) / 6.0);
    } else {
      return Double.MAX_VALUE;
    }
  }

  /**
   * Returns the RMS distance between an arbitrary set of rotamers.
   *
   * @param resrot Residue index-rotamer index pairs.
   * @return RMS distance, or Double.MAX_VALUE if ill-defined.
   */
  public double getRawNBodyDistance(int... resrot) {
    int nRes = resrot.length;
    if (nRes % 2 != 0) {
      throw new IllegalArgumentException(" Must have an even number of arguments; res-rot pairs!");
    }
    nRes /= 2;
    if (nRes < 2) {
      throw new IllegalArgumentException(" Must have >= 4 arguments; at least 2 res-rot pairs!");
    }
    // nCr where n is # of residues, choose pairs.
    int numDists = (int) CombinatoricsUtils.binomialCoefficient(nRes, 2);
    double mult = 1.0 / numDists;
    double totDist2 = 0.0;

    for (int i = 0; i < nRes - 1; i++) {
      int i2 = 2 * i;
      for (int j = i + 1; j < nRes; j++) {
        int j2 = 2 * j;
        double rawDist = checkDistMatrix(resrot[i2], resrot[i2 + 1], resrot[j2], resrot[j2 + 1]);
        if (!Double.isFinite(rawDist) || rawDist == Double.MAX_VALUE) {
          return Double.MAX_VALUE;
        }
        totDist2 += (rawDist * mult);
      }
    }
    return FastMath.sqrt(totDist2);
  }

  /**
   * Checks the distance matrix, finding the shortest distance between the closest rotamers of two
   * residues.
   *
   * @param i Residue i
   * @param ri Rotamer for i
   * @param j Residue j
   * @param rj Rotamer for j
   * @return Shortest distance
   */
  public double getResidueDistance(int i, int ri, int j, int rj) {
    if (i > j) {
      int temp = i;
      i = j;
      j = temp;
      temp = ri;
      ri = rj;
      rj = temp;
    }

    double minDist = Double.MAX_VALUE;
    final int lenri = distanceMatrix[i].length;
    final int lenrj = distanceMatrix[i][ri][j].length;
    for (int roti = 0; roti < lenri; roti++) {
      for (int rotj = 0; rotj < lenrj; rotj++) {
        minDist = Math.min(minDist, checkDistMatrix(i, roti, j, rotj));
      }
    }
    return minDist;
  }

  /**
   * Calculates the minimum distance between two sets of coordinates in a given symmetry operator.
   *
   * @param resi Coordinates of i by [atom][xyz]
   * @param resj Coordinates of j by [atom][xyz]
   * @param symOp Symmetry operator to apply
   * @return Minimum distance
   */
  public double interResidueDistance(double[][] resi, double[][] resj, SymOp symOp) {
    double dist = Double.MAX_VALUE;
    Crystal crystal = molecularAssembly.getCrystal();
    int ni = resi.length;
    for (double[] xi : resi) {
      int nj = resj.length;
      for (double[] xj : resj) {
        if (symOp != null) {
          crystal.applySymOp(xj, xj, symOp);
        }
        // Generally: compare on square-of-distance, and square root only at return.
        // double r = Math.sqrt(crystal.image(xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]));
        double r = crystal.image(xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]);
        if (r < dist) {
          dist = r;
        }
      }
    }
    return sqrt(dist);
  }

  private void distanceMatrix() {

    distanceMatrix = new double[numResidues - 1][][][];
    long numDistances = 0L;
    for (int i = 0; i < (numResidues - 1); i++) {
      Residue residuei = allResiduesArray[i];
      int lengthRi;
      try {
        if (rO.checkIfForced(residuei)) {
          lengthRi = 1;
        } else {
          lengthRi = residuei.getRotamers(library).length;
        }
      } catch (IndexOutOfBoundsException ex) {
        if (useForcedResidues) {
          logger.warning(ex.toString());
        } else {
          rO.logIfMaster(
              format(
                  " Non-forced Residue i %s has null rotamers.",
                  residuei.toFormattedString(false, true)),
              Level.WARNING);
        }
        continue;
      }
      distanceMatrix[i] = new double[lengthRi][][];
      for (int ri = 0; ri < lengthRi; ri++) {
        distanceMatrix[i][ri] = new double[numResidues][];
        for (int j = (i + 1); j < numResidues; j++) {
          Residue residuej = allResiduesArray[j];
          int lengthRj;
          try {
            if (rO.checkIfForced(residuej)) {
              lengthRj = 1;
            } else {
              lengthRj = residuej.getRotamers(library).length;
            }
          } catch (IndexOutOfBoundsException ex) {
            if (useForcedResidues) {
              logger.warning(ex.toString());
            } else {
              rO.logIfMaster(
                  format(
                      " Residue j %s has null rotamers.", residuej.toFormattedString(false, true)));
            }
            continue;
          }
          distanceMatrix[i][ri][j] = new double[lengthRj];
          numDistances += lengthRj;
          if (!lazyMatrix) {
            fill(distanceMatrix[i][ri][j], Double.MAX_VALUE);
          } else {
            fill(distanceMatrix[i][ri][j], -1.0);
          }
        }
      }
    }

    logger.info(format(" Number of pairwise distances: %d", numDistances));

    if (!lazyMatrix) {
      ResidueState[] orig = ResidueState.storeAllCoordinates(allResiduesList);
      int nMultiRes = 0;

      // Build a list that contains one atom from each Residues: CA from
      // amino acids, C1 from nucleic acids, or the first atom otherwise.
      Atom[] atoms = new Atom[numResidues];
      for (int i = 0; i < numResidues; i++) {
        Residue residuei = allResiduesArray[i];
        atoms[i] = residuei.getReferenceAtom();
        if (residuei instanceof MultiResidue) {
          ++nMultiRes;
        }
      }

      /*
       Use of the pre-existing ParallelTeam causes a conflict when
       MultiResidues must re-init the force field. Temporary solution
       for sequence optimization: if > 1 residue optimized, run on only
       one thread.
      */
      int nThreads;
      if (molecularAssembly.getPotentialEnergy().getParallelTeam() != null) {
        nThreads =
            (nMultiRes > 1)
                ? 1
                : molecularAssembly.getPotentialEnergy().getParallelTeam().getThreadCount();
      } else {
        // Suggested: nThreads = (nMultiRes > 1) ? 1 : ParallelTeam.getDefaultThreadCount();
        nThreads = 16;
      }
      ParallelTeam parallelTeam = new ParallelTeam(nThreads);
      Crystal crystal = molecularAssembly.getCrystal();
      int nSymm = crystal.spaceGroup.getNumberOfSymOps();
      logger.info("\n Computing Residue Distance Matrix");

      double nlistCutoff = Math.max(Math.max(distance, twoBodyCutoffDist), threeBodyCutoffDist);

      // Two residues whose c-alphas are separated by 25 angstroms may
      // still interact if they have long side-chains directed at each other.
      double conservativeBuffer = 25.0;
      nlistCutoff += conservativeBuffer;

      /*
       For small crystals, including replicate unit cells in the
       distance matrix is redundant. The cutoff is reduced to the
       interfacial radius.
      */
      if (!crystal.aperiodic()) {
        double sphere =
            min(
                min(crystal.interfacialRadiusA, crystal.interfacialRadiusB),
                crystal.interfacialRadiusC);
        if (nlistCutoff > sphere) {
          nlistCutoff = sphere;
        }
      }

      NeighborList neighborList =
          new NeighborList(null, crystal, atoms, nlistCutoff, 0.0, parallelTeam);

      // Expand coordinates
      double[][] xyz = new double[nSymm][3 * numResidues];
      double[] in = new double[3];
      double[] out = new double[3];
      for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
        SymOp symOp = crystal.spaceGroup.getSymOp(iSymOp);
        for (int i = 0; i < numResidues; i++) {
          int i3 = i * 3;
          int iX = i3;
          int iY = i3 + 1;
          int iZ = i3 + 2;
          Atom atom = atoms[i];
          in[0] = atom.getX();
          in[1] = atom.getY();
          in[2] = atom.getZ();
          crystal.applySymOp(in, out, symOp);
          xyz[iSymOp][iX] = out[0];
          xyz[iSymOp][iY] = out[1];
          xyz[iSymOp][iZ] = out[2];
        }
      }

      // Build the residue neighbor-list.
      int[][][] lists = new int[nSymm][numResidues][];
      boolean[] use = new boolean[numResidues];
      fill(use, true);
      boolean forceRebuild = true;
      boolean printLists = false;
      long neighborTime = -System.nanoTime();
      neighborList.buildList(xyz, lists, use, forceRebuild, printLists);

      neighborTime += System.nanoTime();
      logger.info(
          format(" Built residue neighbor list:           %8.3f sec", neighborTime * 1.0e-9));

      DistanceRegion distanceRegion =
          new DistanceRegion(
              parallelTeam.getThreadCount(),
              numResidues,
              crystal,
              lists,
              neighborList.getPairwiseSchedule());

      long parallelTime = -System.nanoTime();
      try {
        distanceRegion.init(
            this,
            rO,
            molecularAssembly,
            allResiduesArray,
            library,
            algorithmListener,
            distanceMatrix);
        parallelTeam.execute(distanceRegion);
      } catch (Exception e) {
        String message = " Exception compting residue distance matrix.";
        logger.log(Level.SEVERE, message, e);
      }
      parallelTime += System.nanoTime();
      logger.info(
          format(" Pairwise distance matrix:              %8.3f sec\n", parallelTime * 1.0e-9));

      ResidueState.revertAllCoordinates(allResiduesList, orig);

      try {
        parallelTeam.shutdown();
      } catch (Exception ex) {
        logger.warning(
            format(
                " Exception shutting down parallel team for the distance matrix: %s",
                ex.toString()));
      }
    }
  }

  /**
   * Returns the RMS separation distance of three 2-body distances. Defaults to Double.MAX_VALUE
   * when there are pair distances outside cutoffs.
   *
   * @param i Residue i
   * @param ri Rotamer for i
   * @param j Residue j
   * @param rj Rotamer for j
   * @param k Residue k
   * @param rk Rotamer for k
   * @return RMS separation distance
   */
  private double get3BodyDistance(int i, int ri, int j, int rj, int k, int rk) {
    double ij = get2BodyDistance(i, ri, j, rj);
    double ik = get2BodyDistance(i, ri, k, rk);
    double jk = get2BodyDistance(j, rj, k, rk);
    if (areFiniteAndNotMax(ij, ik, jk)) {
      return sqrt((ij * ij + ik * ik + jk * jk) / 3.0);
    } else {
      return Double.MAX_VALUE;
    }
  }

  /**
   * Returns the RMS separation distance of 6 2-body distances. Defaults to Double.MAX_VALUE when
   * there are pair distances outside cutoffs.
   *
   * @param i Residue i
   * @param ri Rotamer for i
   * @param j Residue j
   * @param rj Rotamer for j
   * @param k Residue k
   * @param rk Rotamer for k
   * @param l Residue l
   * @param rl Rotamer for l
   * @return RMS separation distance
   */
  private double get4BodyDistance(int i, int ri, int j, int rj, int k, int rk, int l, int rl) {
    double ij = get2BodyDistance(i, ri, j, rj);
    double ik = get2BodyDistance(i, ri, k, rk);
    double il = get2BodyDistance(i, ri, l, rl);
    double jk = get2BodyDistance(j, rj, k, rk);
    double jl = get2BodyDistance(j, rj, l, rl);
    double kl = get2BodyDistance(k, rk, l, rl);
    if (areFiniteAndNotMax(ij, ik, il, jk, jl, kl)) {
      return sqrt((ij * ij + ik * ik + il * il + jk * jk + jl * jl + kl * kl) / 6.0);
    } else {
      return Double.MAX_VALUE;
    }
  }

  /**
   * Evaluates the pairwise distance between two residues' rotamers under any symmetry operator;
   * does "lazy loading" for the distance matrix.
   *
   * @param i Residue i
   * @param ri Rotamer for i
   * @param j Residue j
   * @param rj Rotamer for j
   * @return Shortest distance
   */
  private double evaluateDistance(int i, int ri, int j, int rj) {
    Residue resi = allResiduesArray[i];
    Rotamer[] rotamersI = resi.getRotamers(library);
    Rotamer roti = rotamersI[ri];
    double[][] xi;
    if (roti.equals(resi.getRotamer())) {
      xi = resi.storeCoordinateArray();
    } else {
      ResidueState origI = resi.storeState();
      RotamerLibrary.applyRotamer(resi, roti);
      xi = resi.storeCoordinateArray();
      resi.revertState(origI);
    }

    Residue resj = allResiduesArray[j];
    Rotamer[] rotamersJ = resj.getRotamers(library);
    Rotamer rotj = rotamersJ[rj];
    double[][] xj;
    if (rotj.equals(resj.getRotamer())) {
      xj = resj.storeCoordinateArray();
    } else {
      ResidueState origJ = resj.storeState();
      RotamerLibrary.applyRotamer(resj, rotj);
      xj = resj.storeCoordinateArray();
      resj.revertState(origJ);
    }

    Crystal crystal = molecularAssembly.getCrystal();
    int nSymm = crystal.spaceGroup.getNumberOfSymOps();
    double minDist = Double.MAX_VALUE;
    for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
      SymOp symOp = crystal.spaceGroup.getSymOp(iSymOp);
      double dist = interResidueDistance(xi, xj, symOp);
      minDist = dist < minDist ? dist : minDist;
    }
    return minDist;
  }
}
