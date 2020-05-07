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
package ffx.potential.nonbonded.implicit;

import static ffx.numerics.atomic.AtomicDoubleArray.atomicDoubleArrayFactory;
import static ffx.numerics.math.DoubleMath.add;
import static ffx.numerics.math.DoubleMath.length2;
import static ffx.numerics.math.DoubleMath.scale;
import static ffx.numerics.math.DoubleMath.sub;
import static java.lang.Double.compare;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import ffx.numerics.atomic.AtomicDoubleArray;
import ffx.numerics.atomic.AtomicDoubleArray.AtomicDoubleArrayImpl;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * GaussVol implements a description molecular volume and surface area described by overlapping
 * Gaussian spheres.
 *
 * <p>Ported from C++ code by Emilio Gallicchio. GaussVol is part of the AGBNP/OpenMM implicit
 * solvent model.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class GaussVol {

  /* -------------------------------------------------------------------------- *
   *                                 GaussVol                                   *
   * -------------------------------------------------------------------------- *
   * This file is part of the AGBNP/OpenMM implicit solvent model software      *
   * implementation funded by the National Science Foundation under grant:      *
   * NSF SI2 1440665  "SI2-SSE: High-Performance Software for Large-Scale       *
   * Modeling of Binding Equilibria"                                            *
   *                                                                            *
   * copyright (c) 2016 Emilio Gallicchio                                       *
   * Authors: Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>                 *
   * Contributors:                                                              *
   *                                                                            *
   *  AGBNP/OpenMM is free software: you can redistribute it and/or modify      *
   *  it under the terms of the GNU Lesser General Public License version 3     *
   *  as published by the Free Software Foundation.                             *
   *                                                                            *
   *  AGBNP/OpenMM is distributed in the hope that it will be useful,           *
   *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
   *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
   *  GNU General Public License for more details.                              *
   *                                                                            *
   *  You should have received a copy of the GNU General Public License         *
   *  along with AGBNP/OpenMM.  If not, see <http://www.gnu.org/licenses/>      *
   *                                                                            *
   * -------------------------------------------------------------------------- */

  private static final Logger logger = Logger.getLogger(GaussVol.class.getName());
  /** Finite-Difference step size to compute surface area. */
  private static final double offset = 0.005;
  /** Conversion factor from a sphere to a Gaussian. */
  private static double KFC = 2.2269859253;

  private static double PFC = 2.5;
  /** Set this to either KFC or PFC. */
  private static double sphereConversion = KFC;
  /** Maximum overlap level. */
  private static int MAX_ORDER = 16;
  /** Volume cutoffs for the switching function (A^3). */
  private static double ANG3 = 1.0;
  /** Volumes smaller than VOLMINB are switched to zero. */
  private static double VOLMINA = 0.01 * ANG3;
  /** Volumes larger than VOLMINB are not switched. */
  private static double VOLMINB = 0.1 * ANG3;
  /** Minimum volume. TODO: Should this be set closer to VOLMINA? */
  private static double MIN_GVOL = Double.MIN_VALUE;
  // private static double MIN_GVOL = 1.0e-12;
  // Output Variables
  private final GaussVolRegion gaussVolRegion;
  private final AtomicDoubleArray3D grad;
  private final SharedDouble totalVolume;
  private final SharedDouble energy;
  private final AtomicDoubleArray gradV;
  private final AtomicDoubleArray freeVolume;
  private final AtomicDoubleArray selfVolume;
  private final ParallelTeam parallelTeam;
  /** Number of atoms. */
  private int nAtoms;
  /**
   * Atomic radii for all atoms. To approximate solvent excluded volume, or solvent accessible
   * surface area, a probe radius can be added to each radius.
   */
  private double[] radii;
  /**
   * All atomic radii with a small offset added, which are used to compute surface area using a
   * finite-difference approach.
   */
  private double[] radiiOffset;
  /** Atomic "self" volumes for all atoms, which are computed from atomic radii. */
  private double[] volumes;
  /** Atomic "self" volumes for all atoms, which are computed from the radii plus offset array. */
  private double[] volumeOffset;
  /** Surface tension -- in our implementation these values are kept at 1.0. */
  private double[] gammas;
  /**
   * Flag to denote if an atom is a hydrogen, which results in it not contributing to the volume or
   * S.A.
   */
  private boolean[] ishydrogen;
  /** The Gaussian Overlap Tree. */
  private GaussianOverlapTree tree;
  /** Maximum depth that the tree reaches */
  private int maximumDepth = 0;
  /** Total number of overlaps in overlap tree */
  private int totalNumberOfOverlaps = 0;
  /** Surface area (Ang^2). */
  private double surfaceArea;
  /** Volume (Ang^3). */
  private double volume;
  /**
   * The volume gradient, which does not include the solvent pressure constant (i.e. this is in
   * Ang^2).
   */
  private double[] volumeGradient;
  /**
   * The surface area gradient, which does not include the surface tension constant (i.e. this is in
   * Ang).
   */
  private double[] surfaceAreaGradient;

  /**
   * Creates/Initializes a GaussVol instance.
   *
   * @param nAtoms The number of atoms.
   * @param radii Atomic radii.
   * @param volumes Atomic volumes.
   * @param gammas Atomic surface tensions.
   * @param ishydrogen True if the atom is a hydrogen.
   * @param parallelTeam ParallelTeam for Parallal jobs.
   */
  public GaussVol(
      int nAtoms,
      double[] radii,
      double[] volumes,
      double[] gammas,
      boolean[] ishydrogen,
      ParallelTeam parallelTeam) {
    tree = new GaussianOverlapTree(nAtoms);
    this.nAtoms = nAtoms;
    this.radii = radii;
    this.volumes = volumes;
    this.gammas = gammas;
    this.ishydrogen = ishydrogen;

    radiiOffset = new double[nAtoms];
    volumeOffset = new double[nAtoms];
    double fourThirdsPI = 4.0 / 3.0 * PI;
    for (int i = 0; i < nAtoms; i++) {
      if (radii[i] != 0.0 && !ishydrogen[i]) {
        radiiOffset[i] = radii[i] + offset;
        volumeOffset[i] = fourThirdsPI * pow(radiiOffset[i], 3);
      } else {
        radiiOffset[i] = 0.0;
        volumeOffset[i] = 0.0;
      }
    }

    this.parallelTeam = parallelTeam;
    int nThreads = 1;
    if (parallelTeam != null) {
      nThreads = parallelTeam.getThreadCount();
      gaussVolRegion = new GaussVolRegion(nThreads);
    } else {
      gaussVolRegion = null;
    }

    totalVolume = new SharedDouble();
    energy = new SharedDouble();
    // Note -- the AtomicDoubleArrayImpl.MULTI does not work yet (use PJ).
    AtomicDoubleArrayImpl atomicDoubleArrayImpl = AtomicDoubleArrayImpl.MULTI;
    grad = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, nThreads);
    gradV = atomicDoubleArrayFactory(atomicDoubleArrayImpl, nThreads, nAtoms);
    freeVolume = atomicDoubleArrayFactory(atomicDoubleArrayImpl, nThreads, nAtoms);
    selfVolume = atomicDoubleArrayFactory(atomicDoubleArrayImpl, nThreads, nAtoms);
  }

  /**
   * Get the radii.
   *
   * @return Returns the radii.
   */
  public double[] getRadii() {
    return radii;
  }

  /**
   * Overlap between two Gaussians represented by a (V,c,a) triplet
   *
   * <p>V: volume of Gaussian c: position of Gaussian a: exponential coefficient
   *
   * <p>g(x) = V (a/pi)^(3/2) exp(-a(x-c)^2)
   *
   * <p>this version is based on V=V(V1,V2,r1,r2,alpha) alpha = (a1 + a2)/(a1 a2)
   *
   * <p>dVdr is (1/r)*(dV12/dr) dVdV is dV12/dV1 dVdalpha is dV12/dalpha d2Vdalphadr is
   * (1/r)*d^2V12/dalpha dr d2VdVdr is (1/r) d^2V12/dV1 dr
   *
   * @param g1 Gaussian 1.
   * @param g2 Gaussian 2.
   * @param g12 Overlap Gaussian.
   * @param dVdr is (1/r)*(dV12/dr)
   * @param dVdV is dV12/dV1
   * @param sfp Derivative of volume.
   * @return Volume.
   */
  private static double overlapGaussianAlpha(
      GaussianVca g1, GaussianVca g2, GaussianVca g12, double[] dVdr, double[] dVdV, double[] sfp) {
    double[] c1 = g1.c;
    double[] c2 = g2.c;
    double[] dist = new double[3];

    sub(c2, c1, dist);
    double d2 = length2(dist);
    double a12 = g1.a + g2.a;
    double deltai = 1.0 / a12;

    // 1/alpha
    double df = (g1.a) * (g2.a) * deltai;
    double ef = exp(-df * d2);
    double gvol = ((g1.v * g2.v) / pow(PI / df, 1.5)) * ef;

    // (1/r)*(dV/dr) w/o switching function
    double dgvol = -2.0 * df * gvol;
    dVdr[0] = dgvol;

    // (1/r)*(dV/dr) w/o switching function
    double dgvolv = g1.v > 0.0 ? gvol / g1.v : 0.0;
    dVdV[0] = dgvolv;

    // Parameters for overlap gaussian.
    // Note that c1 and c2 are Vec3's, and the "*" operator wants the vector first and scalar
    // second:
    // vector2 = vector1 * scalar
    // g12.c = ((c1 * g1.a) + (c2 * g2.a)) * deltai;
    double[] c1a = new double[3];
    double[] c2a = new double[3];
    scale(c1, g1.a * deltai, c1a);
    scale(c2, g2.a * deltai, c2a);
    add(c1a, c2a, g12.c);
    g12.a = a12;
    g12.v = gvol;

    // Switching function
    double[] sp = new double[1];
    double s = switchingFunction(gvol, VOLMINA, VOLMINB, sp);
    sfp[0] = sp[0] * gvol + s;
    return s * gvol;
  }

  /**
   * Overlap volume switching function and 1st derivative.
   *
   * @param gvol Gaussian volume.
   * @param volmina Volume is zero below this limit.
   * @param volminb Volume is not switched above this limit.
   * @param sp Switch first derivative.
   * @return Switch value.
   */
  private static double switchingFunction(
      double gvol, double volmina, double volminb, double[] sp) {
    if (gvol > volminb) {
      sp[0] = 0.0;
      return 1.0;
    } else if (gvol < volmina) {
      sp[0] = 0.0;
      return 0.0;
    }
    double swd = 1.0 / (volminb - volmina);
    double swu = (gvol - volmina) * swd;
    double swu2 = swu * swu;
    double swu3 = swu * swu2;
    sp[0] = swd * 30.0 * swu2 * (1.0 - 2.0 * swu + swu2);
    return swu3 * (10.0 - 15.0 * swu + 6.0 * swu2);
  }

  /**
   * Compute molecular volume and surface area.
   *
   * @param positions Atomic positions to use.
   * @return The volume.
   */
  public double computeVolumeAndSA(double[][] positions) {

    if (parallelTeam == null || parallelTeam.getThreadCount() == 1) {
      // Update the overlap tree.
      computeTree(positions);

      // Compute the volume.
      computeVolume(totalVolume, energy, grad, gradV, freeVolume, selfVolume);
    } else {
      // Execute in parallel.
      try {
        gaussVolRegion.init(0, positions);
        parallelTeam.execute(gaussVolRegion);
      } catch (Exception e) {
        logger.severe(" Exception evaluating GaussVol " + e.toString());
      }
    }

    if (volumeGradient == null || volumeGradient.length != nAtoms * 3) {
      volumeGradient = new double[nAtoms * 3];
      surfaceAreaGradient = new double[nAtoms * 3];
    }

    double selfVolumeSum = 0;
    for (int i = 0; i < nAtoms; i++) {
      selfVolumeSum += selfVolume.get(i);
      double x = grad.getX(i);
      double y = grad.getY(i);
      double z = grad.getZ(i);
      if (ishydrogen[i]) {
        x = 0.0;
        y = 0.0;
        z = 0.0;
      }
      // Volume based gradient
      int index = i * 3;
      volumeGradient[index++] = x;
      volumeGradient[index++] = y;
      volumeGradient[index] = z;

      // Surface Area based gradient
      index = i * 3;
      surfaceAreaGradient[index++] = x;
      surfaceAreaGradient[index++] = y;
      surfaceAreaGradient[index] = z;
    }

    // Save the total molecular volume.
    volume = totalVolume.get();

    // Save a reference to non-offset radii and volumes.
    double[] radiiBak = this.radii;
    double[] volumeBak = this.volumes;

    // Load the offset radii and volumes.
    this.radii = radiiOffset;
    this.volumes = volumeOffset;

    // Run Volume calculation on radii that are slightly offset in order to do finite difference to
    // get back surface area
    if (parallelTeam == null || parallelTeam.getThreadCount() == 1) {
      // Rescan the overlap tree.
      rescanTreeVolumes(positions);
      // Compute the volume.
      computeVolume(totalVolume, energy, grad, gradV, freeVolume, selfVolume);
    } else {
      // Execute in parallel.
      try {
        gaussVolRegion.init(2, positions);
        parallelTeam.execute(gaussVolRegion);
      } catch (Exception e) {
        logger.severe(" Exception evaluating GaussVol " + e.toString());
      }
    }

    // Set the radii and volumes to their non-offset values.
    this.radii = radiiBak;
    this.volumes = volumeBak;

    double selfVolumeOffsetSum = 0;
    for (int i = 0; i < nAtoms; i++) {
      selfVolumeOffsetSum += selfVolume.get(i);
      double x = grad.getX(i);
      double y = grad.getY(i);
      double z = grad.getZ(i);
      if (ishydrogen[i]) {
        x = 0.0;
        y = 0.0;
        z = 0.0;
      }
      // Surface Area based gradient
      int index = i * 3;
      surfaceAreaGradient[index] = (x - surfaceAreaGradient[index++]) / offset;
      surfaceAreaGradient[index] = (y - surfaceAreaGradient[index++]) / offset;
      surfaceAreaGradient[index] = (z - surfaceAreaGradient[index]) / offset;
    }

    // Calculate the surface area.
    surfaceArea = (selfVolumeOffsetSum - selfVolumeSum) / offset;

    return volume;
  }

  /**
   * Returns the maximum depth of the overlap tree
   *
   * @return maximumDepth
   */
  public int getMaximumDepth() {
    return maximumDepth;
  }

  /**
   * Return Surface Area (A^2).
   *
   * @return Surface Area (A^2).
   */
  public double getSurfaceArea() {
    return surfaceArea;
  }

  public double[] getSurfaceAreaGradient() {
    return surfaceAreaGradient;
  }

  /**
   * Return the total number of overlaps in the tree
   *
   * @return totalNumberOfOverlaps
   */
  public int getTotalNumberOfOverlaps() {
    return totalNumberOfOverlaps;
  }

  /**
   * Return Volume (A^3).
   *
   * @return Volume (A^3).
   */
  public double getVolume() {
    return volume;
  }

  public double[] getVolumeGradient() {
    return volumeGradient;
  }

  /**
   * Set gamma values.
   *
   * @param gammas Gamma values (kcal/mol/A^2).
   * @throws Exception If the array length does not match the number of atoms.
   */
  public void setGammas(double[] gammas) throws Exception {
    if (nAtoms == gammas.length) {
      this.gammas = gammas;
    } else {
      throw new Exception(" setGammas: number of atoms does not match");
    }
  }

  /**
   * Set the isHydrogen flag.
   *
   * @param isHydrogen Per atom flag to indicate if its a hydrogen.
   * @throws Exception If the array length does not match the number of atoms.
   */
  public void setIsHydrogen(boolean[] isHydrogen) throws Exception {
    if (nAtoms == isHydrogen.length) {
      this.ishydrogen = isHydrogen;
    } else {
      throw new Exception(" setIsHydrogen: number of atoms does not match");
    }
  }

  /**
   * Set radii and volumes.
   *
   * @param radii Atomic radii (Angstroms).
   * @param volumes Atomic volumes (Angstroms^3).
   * @throws Exception If the array length does not match the number of atoms.
   */
  public void setRadiiAndVolumes(double[] radii, double[] volumes) throws Exception {
    if (nAtoms == radii.length) {
      this.radii = radii;
      radiiOffset = new double[nAtoms];
      volumeOffset = new double[nAtoms];
      double fourThirdsPI = 4.0 / 3.0 * PI;
      for (int i = 0; i < nAtoms; i++) {
        radiiOffset[i] = radii[i] + offset;
        volumeOffset[i] = fourThirdsPI * pow(radiiOffset[i], 3);
      }
    } else {
      throw new Exception(" setRadii: number of atoms does not match");
    }

    if (nAtoms == volumes.length) {
      this.volumes = volumes;
    } else {
      throw new Exception(" setVolumes: number of atoms does not match");
    }
  }

  /**
   * Constructs the tree.
   *
   * @param positions Current atomic positions.
   */
  private void computeTree(double[][] positions) {
    tree.computeOverlapTreeR(positions, radii, volumes, gammas, ishydrogen);
  }

  /**
   * Returns GaussVol volume energy function and forces. Also returns gradients with respect to
   * atomic volumes and atomic free-volumes and self-volumes.
   *
   * @param totalVolume Total volume.
   * @param totalEnergy Total energy.
   * @param grad Atomic gradient.
   * @param gradV Volume gradient.
   * @param freeVolume Free volume.
   * @param selfVolume Self volume.
   */
  private void computeVolume(
      SharedDouble totalVolume,
      SharedDouble totalEnergy,
      AtomicDoubleArray3D grad,
      AtomicDoubleArray gradV,
      AtomicDoubleArray freeVolume,
      AtomicDoubleArray selfVolume) {
    // Initialize output variables.
    totalVolume.set(0.0);
    totalEnergy.set(0.0);
    grad.reset(0, 0, nAtoms - 1);
    gradV.reset(0, 0, nAtoms - 1);
    freeVolume.reset(0, 0, nAtoms - 1);
    selfVolume.reset(0, 0, nAtoms - 1);

    // Compute the volume.
    tree.computeVolume2R(totalVolume, totalEnergy, grad, gradV, freeVolume, selfVolume);

    // Reduce output variables.
    grad.reduce(0, nAtoms - 1);
    gradV.reduce(0, nAtoms - 1);
    freeVolume.reduce(0, nAtoms - 1);
    selfVolume.reduce(0, nAtoms - 1);

    for (int i = 0; i < nAtoms; ++i) {
      if (volumes[i] > 0) {
        gradV.set(0, i, gradV.get(i) / volumes[i]);
      }
    }
  }

  /**
   * Rescan the tree after resetting gammas, radii and volumes.
   *
   * @param positions Atomic coordinates.
   */
  private void rescanTreeVolumes(double[][] positions) {
    tree.rescanTreeV(positions, radii, volumes, gammas, ishydrogen);
  }

  /** Rescan the tree resetting gammas only with current values. */
  private void rescanTreeGammas() {
    tree.rescanTreeG(gammas);
  }

  /**
   * Returns number of overlaps for each atom.
   *
   * @param nov Number of overlaps.
   */
  private void getStats(int[] nov) {
    if (nov.length != nAtoms) return;

    for (int i = 0; i < nAtoms; i++) nov[i] = 0;
    for (int atom = 0; atom < nAtoms; atom++) {
      int slot = atom + 1;
      nov[atom] = tree.nChildrenUnderSlotR(slot);
    }
  }

  /** Print the tree. */
  private void printTree() {
    tree.printTree();
  }

  /** 3D Gaussian, V,c,a representation. */
  private static class GaussianVca {
    // Gaussian volume
    public double v;
    // Gaussian exponent
    public double a;
    // Center
    public double[] c = new double[3];
  }

  /**
   * Overlap between two Gaussians represented by a (V,c,a) triplet V: volume of Gaussian c:
   * position of Gaussian a: exponential coefficient
   *
   * <p>g(x) = V (a/pi)^(3/2) exp(-a(x-c)^2)
   *
   * <p>This version is based on V=V(V1,V2,r1,r2,alpha) alpha = (a1 + a2)/(a1 a2) dVdr is
   * (1/r)*(dV12/dr) dVdV is dV12/dV1 dVdalpha is dV12/dalpha d2Vdalphadr is (1/r)*d^2V12/dalpha dr
   * d2VdVdr is (1/r) d^2V12/dV1 dr
   */
  private static class GaussianOverlap implements Comparable<GaussianOverlap> {
    /** level (0=root, 1=atoms, 2=2-body, 3=3-body, etc.) */
    public int level;
    /** Gaussian representing overlap */
    public GaussianVca g;
    /** Volume of overlap (also stores Psi1..i in GPU version) */
    public double volume;
    /** The atomic index of the last atom of the overlap list (i, j, k, ..., atom) */
    public int atom;
    /** Derivative wrt volume of first atom (also stores F1..i in GPU version) */
    double dvv1;
    /** Derivative wrt position of first atom (also stores P1..i in GPU version) */
    double[] dv1;
    /** Sum gammai for this overlap (surface tension parameter) */
    double gamma1i;
    /** Self volume accumulator (also stores Psi'1..i in GPU version) */
    double selfVolume;
    /** Switching function derivatives */
    double sfp;
    /** = (Parent, atom) index in tree list of parent overlap */
    int parentIndex;
    /** Start index in tree array of children */
    int childrenStartIndex;
    /** Number of children. */
    int childrenCount;

    public GaussianOverlap() {
      g = new GaussianVca();
      dv1 = new double[3];
    }

    public GaussianOverlap(GaussianVca g, double volume, double selfVolume, int atom) {
      this.g = g;
      this.volume = volume;
      this.selfVolume = selfVolume;
      this.atom = atom;

      dv1 = new double[3];
    }

    /**
     * Order by volume, larger first.
     *
     * @param o GaussianOverlap to compare to.
     * @return Compare by volume.
     */
    @Override
    public int compareTo(GaussianOverlap o) {
      return compare(volume, o.volume);
    }

    /** Print overlaps. */
    public void printOverlap() {
      logger.info(
          format(
              " Gaussian Overlap %d: Atom: %d, Parent: %d, ChildrenStartIndex: %d, ChildrenCount: %d,"
                  + "Volume: %6.3f, Gamma: %6.3f, Gauss.a: %6.3f, Gauss.v: %6.3f, Gauss.center (%6.3f,%6.3f,%6.3f),"
                  + "dedx: %6.3f, dedy: %6.3f, dedz: %6.3f, sfp: %6.3f",
              level,
              atom,
              parentIndex,
              childrenStartIndex,
              childrenCount,
              volume,
              gamma1i,
              g.a,
              g.v,
              g.c[0],
              g.c[1],
              g.c[2],
              dv1[0],
              dv1[1],
              dv1[2],
              sfp));
    }
  }

  /** Gaussian Overlap Tree. */
  private class GaussianOverlapTree {

    /** Number of atoms. */
    int nAtoms;
    /** The root is at index 0. Atoms are from 1 .. nAtoms. */
    List<GaussianOverlap> overlaps;

    /**
     * GaussianOverlapTree constructor.
     *
     * @param nAtoms Number of atoms.
     */
    GaussianOverlapTree(int nAtoms) {
      this.nAtoms = nAtoms;
      overlaps = Collections.synchronizedList(new ArrayList<>(nAtoms + 1));
    }

    /**
     * Initialize the Overlap Tree.
     *
     * @param pos Atomic positions.
     * @param radii Atomic radii.
     * @param volumes Atomic volumes.
     * @param gammas Atomic surface tensions.
     * @param ishydrogen True if the atom is a hydrogen.
     */
    void initOverlapTree(
        double[][] pos, double[] radii, double[] volumes, double[] gammas, boolean[] ishydrogen) {

      // Reset tree
      overlaps = Collections.synchronizedList(new ArrayList<>(nAtoms + 1));

      // Slot 0 contains the master tree information, children = all of the atoms.
      GaussianOverlap overlap = new GaussianOverlap();
      overlap.level = 0;
      overlap.volume = 0;
      overlap.dvv1 = 0.;
      overlap.selfVolume = 0;
      overlap.sfp = 1.;
      overlap.gamma1i = 0.;
      overlap.parentIndex = -1;
      overlap.atom = -1;
      overlap.childrenStartIndex = 1;
      overlap.childrenCount = nAtoms;
      overlaps.add(0, overlap);

      // List of atoms start at slot 1.
      for (int iat = 0; iat < nAtoms; iat++) {
        overlap = new GaussianOverlap();
        double a = sphereConversion / (radii[iat] * radii[iat]);
        double vol = ishydrogen[iat] ? 0 : volumes[iat];
        overlap.level = 1;
        overlap.g.v = vol;
        overlap.g.a = a;
        overlap.g.c = pos[iat];
        overlap.volume = vol;
        overlap.dvv1 = 1; // dVi/dVi
        overlap.selfVolume = 0;
        overlap.sfp = 1;
        overlap.gamma1i = gammas[iat]; // gamma[iat] / SA_DR;
        overlap.parentIndex = 0;
        overlap.atom = iat;
        overlap.childrenStartIndex = -1;
        overlap.childrenCount = -1;
        overlaps.add(iat + 1, overlap);
      }
    }

    /**
     * Add Children.
     *
     * @param parentIndex Parent index.
     * @param childrenOverlaps Children overlaps.
     * @return Index of the first added child.
     */
    int addChildren(int parentIndex, List<GaussianOverlap> childrenOverlaps) {

      // Adds children starting at the last slot
      int startIndex = overlaps.size();

      int noverlaps = childrenOverlaps.size();

      // Retrieves address of root overlap
      GaussianOverlap root = overlaps.get(parentIndex);

      // Registers list of children
      root.childrenStartIndex = startIndex;
      root.childrenCount = noverlaps;

      // Sort neighbors by overlap volume.
      Collections.sort(childrenOverlaps);

      int rootLevel = root.level;
      int nextLevel = rootLevel + 1;

      // Now copies the children overlaps from temp buffer.
      for (GaussianOverlap child : childrenOverlaps) {
        child.level = nextLevel;

        // Connect overlap to parent
        child.parentIndex = parentIndex;

        // Reset its children indexes
        child.childrenStartIndex = -1;
        child.childrenCount = -1;

        // Add to tree.
        overlaps.add(child);
      }

      return startIndex;
    }

    /**
     * Scans the siblings of overlap identified by "rootIndex" to create children overlaps, returns
     * them into the "childrenOverlaps" buffer: (root) + (atom) -> (root, atom)
     *
     * @param rootIndex Root index.
     * @param childrenOverlaps Children overlaps.
     */
    void computeChildren(int rootIndex, List<GaussianOverlap> childrenOverlaps) {
      int parentIndex;
      int siblingStart, siblingCount;

      // Reset output buffer
      childrenOverlaps.clear();

      // Retrieves overlap.
      GaussianOverlap root = overlaps.get(rootIndex);

      // Retrieves parent overlap.
      parentIndex = root.parentIndex;

      // Master root? can't do computeChildren() on master root
      if (parentIndex < 0) {
        throw new IllegalArgumentException(" Cannot compute children of master node!");
      }

      if (root.level >= MAX_ORDER) {
        return; // Ignore overlaps above a certain order to cap computational cost.
      }

      GaussianOverlap parent = overlaps.get(parentIndex);

      // Retrieves start index and count of siblings. Includes both younger and older siblings
      siblingStart = parent.childrenStartIndex;
      siblingCount = parent.childrenCount;

      // Parent is not initialized?
      if (siblingStart < 0 || siblingCount < 0) {
        throw new IllegalArgumentException(
            format(" Parent %s of overlap %s has no sibilings.", parent, root));
      }

      // This overlap somehow is not the child of registered parent.
      if (rootIndex < siblingStart && rootIndex > siblingStart + siblingCount - 1) {
        throw new IllegalArgumentException(
            format(" Node %s is somehow not the child of its parent %s", root, parent));
      }

      // Now loops over "younger" siblings (i<j loop) to compute new overlaps.
      // Loop starts at the first younger sibling, and runs to the end of all siblings.
      for (int slotj = rootIndex + 1; slotj < siblingStart + siblingCount; slotj++) {

        GaussianVca g12 = new GaussianVca();

        GaussianOverlap sibling = overlaps.get(slotj);
        double gvol;
        double[] dVdr = new double[1];
        double[] dVdV = new double[1];
        double[] sfp = new double[1];

        // Atomic gaussian of last atom of sibling.
        int atom2 = sibling.atom;
        GaussianVca g1 = root.g;

        // Atoms are stored in the tree at indexes 1...N
        GaussianVca g2 = overlaps.get(atom2 + 1).g;
        gvol = overlapGaussianAlpha(g1, g2, g12, dVdr, dVdV, sfp);

        /*
         Create child if overlap volume is above a threshold.
         Due to Gaussians having infinite support, volume is never zero.
        */
        if (gvol > MIN_GVOL) {
          GaussianOverlap ov = new GaussianOverlap(g12, gvol, 0.0, atom2);
          // dv1 is the gradient of V(123..)n with respect to the position of 1
          // ov.dv1 = ( g2.c - g1.c ) * (-dVdr);
          sub(g2.c, g1.c, ov.dv1);
          scale(ov.dv1, -dVdr[0], ov.dv1);

          // dvv1 is the derivative of V(123...)n with respect to V(123...)
          ov.dvv1 = dVdV[0];
          ov.sfp = sfp[0];
          ov.gamma1i = root.gamma1i + overlaps.get(atom2 + 1).gamma1i;
          childrenOverlaps.add(ov);
        }
      }
    }

    /**
     * Grow the tree with more children starting at the given root slot (recursive).
     *
     * @param root The root index.
     */
    private void computeAndAddChildrenR(int root) {
      List<GaussianOverlap> childrenOverlaps = Collections.synchronizedList(new ArrayList<>());
      computeChildren(root, childrenOverlaps);
      int nOverlaps = childrenOverlaps.size();
      if (nOverlaps > 0) {
        int startSlot = addChildren(root, childrenOverlaps);
        for (int ichild = startSlot; ichild < startSlot + nOverlaps; ichild++) {
          computeAndAddChildrenR(ichild);
          totalNumberOfOverlaps++;
        }
      }
    }

    /**
     * Compute the overlap tree.
     *
     * @param pos Atomic positions.
     * @param radii Atomic radii.
     * @param volumes Atomic volumes.
     * @param gammas Atomic surface tensions.
     * @param ishydrogen True if the atom is a hydrogen.
     */
    void computeOverlapTreeR(
        double[][] pos, double[] radii, double[] volumes, double[] gammas, boolean[] ishydrogen) {
      initOverlapTree(pos, radii, volumes, gammas, ishydrogen);
      for (int slot = 1; slot <= nAtoms; slot++) {
        computeAndAddChildrenR(slot);
        totalNumberOfOverlaps++;
      }
    }

    /**
     * Compute volumes, energy of the overlap at slot and calls itself recursively to get the
     * volumes of the children.
     *
     * @param slot Slot to begin from.
     * @param psi1i Subtree accumulator for free volume.
     * @param f1i Subtree accumulator for free volume.
     * @param p1i Subtree accumulator for free volume.
     * @param psip1i Subtree accumulator for self volume.
     * @param fp1i Subtree accumulators for self volume.
     * @param pp1i Subtree accumulators for self volume.
     * @param energy1i Subtree accumulator for volume-based energy.
     * @param fenergy1i Subtree accumulator for volume-based energy.
     * @param penergy1i subtree accumulator for volume-based energy.
     * @param threadID Thread for accumulation.
     * @param dr Gradient of volume-based energy wrt to atomic positions.
     * @param dv Gradient of volume-based energy wrt to atomic volumes.
     * @param freeVolume Atomic free volumes.
     * @param selfVolume Atomic self volumes.
     */
    void computeVolumeUnderSlot2R(
        int slot,
        double[] psi1i,
        double[] f1i,
        double[] p1i,
        double[] psip1i,
        double[] fp1i,
        double[] pp1i,
        double[] energy1i,
        double[] fenergy1i,
        double[] penergy1i,
        int threadID,
        AtomicDoubleArray3D dr,
        AtomicDoubleArray dv,
        AtomicDoubleArray freeVolume,
        AtomicDoubleArray selfVolume) {

      GaussianOverlap ov = overlaps.get(slot);
      // Keep track of overlap depth for each overlap.
      // If a new depth is greater than previous greatest, save depth in maximumDepth
      if (ov.level >= maximumDepth) {
        maximumDepth = ov.level;
      }

      // Whether to add volumes (e.g. selfs, 3-body overlaps) or subtract them (e.g. 2-body
      // overlaps, 4-body overlaps)
      double cf = ov.level % 2 == 0 ? -1.0 : 1.0;
      // Overall volume is increased/decreased by the full volume.
      double volcoeff = ov.level > 0 ? cf : 0;
      // Atomic contributions to overlap volume are evenly distributed.
      double volcoeffp = ov.level > 0 ? volcoeff / (double) ov.level : 0;

      int atom = ov.atom;
      double ai = overlaps.get(atom + 1).g.a;
      double a1i = ov.g.a;
      double a1 = a1i - ai;

      // For free volumes
      psi1i[0] = volcoeff * ov.volume;
      f1i[0] = volcoeff * ov.sfp;

      // For self volumes
      psip1i[0] = volcoeffp * ov.volume;
      fp1i[0] = volcoeffp * ov.sfp;

      // EV energy
      energy1i[0] = volcoeffp * ov.gamma1i * ov.volume;
      fenergy1i[0] = volcoeffp * ov.sfp * ov.gamma1i;

      // Loop over children.
      if (ov.childrenStartIndex >= 0) {
        for (int sloti = ov.childrenStartIndex;
            sloti < ov.childrenStartIndex + ov.childrenCount;
            sloti++) {
          double[] psi1it = new double[1];
          double[] f1it = new double[1];
          double[] p1it = new double[3];
          double[] psip1it = new double[1];
          double[] fp1it = new double[1];
          double[] pp1it = new double[3];
          double[] energy1it = new double[1];
          double[] fenergy1it = new double[1];
          double[] penergy1it = new double[3];
          computeVolumeUnderSlot2R(
              sloti,
              psi1it,
              f1it,
              p1it,
              psip1it,
              fp1it,
              pp1it,
              energy1it,
              fenergy1it,
              penergy1it,
              threadID,
              dr,
              dv,
              freeVolume,
              selfVolume);
          psi1i[0] += psi1it[0];
          f1i[0] += f1it[0];
          add(p1i, p1it, p1i);

          psip1i[0] += psip1it[0];
          fp1i[0] += fp1it[0];
          add(pp1i, pp1it, pp1i);

          energy1i[0] += energy1it[0];
          fenergy1i[0] += fenergy1it[0];
          add(penergy1i, penergy1it, penergy1i);
        }
      }

      // Skip this for the Root Level.
      if (ov.level > 0) {
        // Contributions to free and self volume of last atom
        freeVolume.add(threadID, atom, psi1i[0]);
        selfVolume.add(threadID, atom, psip1i[0]);

        // Contributions to energy gradients
        double c2 = ai / a1i;

        // dr[atom] += (-ov.dv1) * fenergy1i + penergy1i * c2;
        double[] work1 = new double[3];
        double[] work2 = new double[3];
        double[] work3 = new double[3];
        scale(penergy1i, c2, work1);
        scale(ov.dv1, -fenergy1i[0], work2);
        add(work1, work2, work3);
        dr.add(threadID, atom, work3[0], work3[1], work3[2]);

        // ov.g.v is the unswitched volume
        dv.add(threadID, atom, ov.g.v * fenergy1i[0]);

        // Update subtree P1..i's for parent
        c2 = a1 / a1i;

        // p1i = (ov.dv1) * f1i + p1i * c2;
        scale(ov.dv1, f1i[0], work1);
        scale(pp1i, c2, work2);
        add(work1, work2, p1i);

        // pp1i = (ov.dv1) * fp1i + pp1i * c2;
        scale(ov.dv1, fp1i[0], work1);
        scale(pp1i, c2, work2);
        add(work1, work2, pp1i);

        // penergy1i = (ov.dv1) * fenergy1i + penergy1i * c2;
        scale(ov.dv1, fenergy1i[0], work1);
        scale(penergy1i, c2, work2);
        add(work1, work2, penergy1i);

        // Update subtree F1..i's for parent
        f1i[0] = ov.dvv1 * f1i[0];
        fp1i[0] = ov.dvv1 * fp1i[0];
        fenergy1i[0] = ov.dvv1 * fenergy1i[0];
      }
    }

    /**
     * Recursively traverses tree and computes volumes, etc.
     *
     * @param volume Volume.
     * @param energy Energy.
     * @param dr Coordinate derivatives.
     * @param dv Volume derivatives.
     * @param freeVolume Free volume.
     * @param selfvolume Self volume.
     */
    void computeVolume2R(
        SharedDouble volume,
        SharedDouble energy,
        AtomicDoubleArray3D dr,
        AtomicDoubleArray dv,
        AtomicDoubleArray freeVolume,
        AtomicDoubleArray selfvolume) {

      double[] psi1i = new double[1];
      double[] f1i = new double[1];
      double[] p1i = new double[3];
      double[] psip1i = new double[1];
      double[] fp1i = new double[1];
      double[] pp1i = new double[3];
      double[] energy1i = new double[1];
      double[] fenergy1i = new double[1];
      double[] penergy1i = new double[3];
      // Only one thread for serial computation.
      int threadID = 0;

      computeVolumeUnderSlot2R(
          0,
          psi1i,
          f1i,
          p1i,
          psip1i,
          fp1i,
          pp1i,
          energy1i,
          fenergy1i,
          penergy1i,
          threadID,
          dr,
          dv,
          freeVolume,
          selfvolume);

      volume.addAndGet(psi1i[0]);
      energy.addAndGet(energy1i[0]);
    }

    /**
     * Rescan the tree to recompute the volumes. It does not modify the tree.
     *
     * @param pos Atomic positions.
     * @param radii Atomic radii.
     * @param volumes Atomic volumes.
     * @param gammas Atomic surface tensions.
     * @param ishydrogen True if the atom is a hydrogen.
     */
    void rescanTreeV(
        double[][] pos, double[] radii, double[] volumes, double[] gammas, boolean[] ishydrogen) {
      initRescanTreeV(pos, radii, volumes, gammas, ishydrogen);
      rescanR(0);
    }

    /**
     * Rescan the sub-tree to recompute the volumes. It does not modify the tree.
     *
     * @param slot The slot to begin from.
     */
    void rescanR(int slot) {
      int parentIndex;

      // This overlap.
      GaussianOverlap ov = overlaps.get(slot);

      // Recompute its own overlap by merging parent and last atom.
      parentIndex = ov.parentIndex;
      if (parentIndex > 0) {
        GaussianVca g12 = new GaussianVca();
        double[] dVdr = new double[1];
        double[] dVdV = new double[1];
        double[] sfp = new double[1];

        int atom = ov.atom;
        GaussianOverlap parent = overlaps.get(parentIndex);
        GaussianVca g1 = parent.g;

        // Atoms are stored in the tree at indexes 1...N
        GaussianVca g2 = overlaps.get(atom + 1).g;
        double gvol = overlapGaussianAlpha(g1, g2, g12, dVdr, dVdV, sfp);
        ov.g = g12;
        ov.volume = gvol;

        // dv1 is the gradient of V(123..)n with respect to the position of 1
        // ov.dv1 = ( g2.c - g1.c ) * (-dVdr);
        sub(g2.c, g1.c, ov.dv1);
        scale(ov.dv1, -dVdr[0], ov.dv1);

        // dvv1 is the derivative of V(123...)n with respect to V(123...)
        ov.dvv1 = dVdV[0];
        ov.sfp = sfp[0];
        ov.gamma1i = parent.gamma1i + overlaps.get(atom + 1).gamma1i;
      }

      // Calls itself recursively on the children.
      for (int slotChild = ov.childrenStartIndex;
          slotChild < ov.childrenStartIndex + ov.childrenCount;
          slotChild++) {
        rescanR(slotChild);
      }
    }

    /**
     * Init rescan of the tree to recompute the volumes. It does not modify the tree.
     *
     * @param pos Atomic coodinates.
     * @param radii Atomic radii.
     * @param volumes Atomic volumes.
     * @param gammas Atomic surface tensions.
     * @param ishydrogen True if the atom is a hydrogen.
     */
    void initRescanTreeV(
        double[][] pos, double[] radii, double[] volumes, double[] gammas, boolean[] ishydrogen) {
      int slot = 0;
      GaussianOverlap ov = overlaps.get(slot);
      ov.level = 0;
      ov.volume = 0.0;
      ov.dv1 = new double[3];
      ov.dvv1 = 0.0;
      ov.selfVolume = 0.0;
      ov.sfp = 1.0;
      ov.gamma1i = 0.0;

      slot = 1;
      for (int iat = 0; iat < nAtoms; iat++, slot++) {
        double a = sphereConversion / (radii[iat] * radii[iat]);
        double vol = ishydrogen[iat] ? 0.0 : volumes[iat];
        ov = overlaps.get(slot);
        ov.level = 1;
        ov.g.v = vol;
        ov.g.a = a;
        ov.g.c = pos[iat];
        ov.volume = vol;
        ov.dv1 = new double[3];
        ov.dvv1 = 1.0; // dVi/dVi
        ov.selfVolume = 0.0;
        ov.sfp = 1.0;
        ov.gamma1i = gammas[iat]; // gamma[iat]/SA_DR
      }
    }

    /**
     * Rescan the sub-tree to recompute the gammas. It does not modify the volumes nor the tree.
     *
     * @param slot Slot to begin from.
     */
    void rescanGammaR(int slot) {

      // This overlap.
      GaussianOverlap go = overlaps.get(slot);

      // Recompute its own overlap by merging parent and last atom.
      int parentIndex = go.parentIndex;
      if (parentIndex > 0) {
        int atom = go.atom;
        GaussianOverlap parent = overlaps.get(parentIndex);
        go.gamma1i = parent.gamma1i + overlaps.get(atom + 1).gamma1i;
      }

      // Calls itself recursively on the children.
      for (int slotChild = go.childrenStartIndex;
          slotChild < go.childrenStartIndex + go.childrenCount;
          slotChild++) {
        rescanGammaR(slotChild);
      }
    }

    /**
     * Rescan the tree to recompute the gammas only. It does not modify volumes and the tree.
     *
     * @param gammas Gamma values.
     */
    void rescanTreeG(double[] gammas) {

      int slot = 0;
      GaussianOverlap ov = overlaps.get(slot);
      ov.gamma1i = 0.;

      slot = 1;
      for (int iat = 0; iat < nAtoms; iat++, slot++) {
        ov = overlaps.get(slot);
        ov.gamma1i = gammas[iat];
      }

      rescanGammaR(0);
    }

    /** Print the contents of the tree. */
    void printTree() {
      // logger.info("slot level LastAtom parent ChStart ChCount SelfV V gamma a x y z dedx dedy
      // dedz sfp");
      for (int i = 1; i <= nAtoms; i++) {
        printTreeR(i);
      }
    }

    /**
     * Print the contents of the tree (recursive).
     *
     * @param slot Slot to begin from.
     */
    void printTreeR(int slot) {
      GaussianOverlap ov = overlaps.get(slot);
      logger.info(format("tg:      %d ", slot));
      ov.printOverlap();
      for (int i = ov.childrenStartIndex; i < ov.childrenStartIndex + ov.childrenCount; i++) {
        printTreeR(i);
      }
    }

    /**
     * Counts number of overlaps under the one given.
     *
     * @param slot Slot to begin from.
     * @return Number of overlaps.
     */
    int nChildrenUnderSlotR(int slot) {
      int n = 0;
      if (overlaps.get(slot).childrenCount > 0) {
        n += overlaps.get(slot).childrenCount;
        // now calls itself on the children
        for (int i = 0; i < overlaps.get(slot).childrenCount; i++) {
          n += nChildrenUnderSlotR(overlaps.get(slot).childrenStartIndex + i);
        }
      }
      return n;
    }
  }

  /** Use instances of the GaussianOverlapTree to compute the GaussVol volume in parallel. */
  private class GaussVolRegion extends ParallelRegion {

    private GaussianOverlapTree[] localTree;
    private ComputeTreeLoop[] computeTreeLoops;
    private ComputeVolumeLoop[] computeVolumeLoops;
    private RescanTreeLoop[] rescanTreeLoops;
    private ReductionLoop[] reductionLoops;
    private int mode = 0;
    private double[][] coordinates = null;

    /**
     * GaussVolRegion constructor.
     *
     * @param nThreads Number of threads.
     */
    public GaussVolRegion(int nThreads) {
      localTree = new GaussianOverlapTree[nThreads];
      computeTreeLoops = new ComputeTreeLoop[nThreads];
      computeVolumeLoops = new ComputeVolumeLoop[nThreads];
      rescanTreeLoops = new RescanTreeLoop[nThreads];
      reductionLoops = new ReductionLoop[nThreads];
    }

    public void init(int mode, double[][] coordinates) {
      this.mode = mode;
      this.coordinates = coordinates;

      // Initialize output variables.
      totalVolume.set(0.0);
      energy.set(0.0);
      grad.reset(parallelTeam, 0, nAtoms - 1);
      gradV.reset(parallelTeam, 0, nAtoms - 1);
      freeVolume.reset(parallelTeam, 0, nAtoms - 1);
      selfVolume.reset(parallelTeam, 0, nAtoms - 1);
    }

    @Override
    public void run() throws Exception {
      int threadIndex = getThreadIndex();
      if (computeTreeLoops[threadIndex] == null) {
        localTree[threadIndex] = new GaussianOverlapTree(nAtoms);
        computeTreeLoops[threadIndex] = new ComputeTreeLoop();
        computeVolumeLoops[threadIndex] = new ComputeVolumeLoop();
        rescanTreeLoops[threadIndex] = new RescanTreeLoop();
        reductionLoops[threadIndex] = new ReductionLoop();
      }
      try {
        if (mode == 0) {
          execute(1, nAtoms, computeTreeLoops[threadIndex]);
          execute(1, nAtoms, computeVolumeLoops[threadIndex]);
          execute(0, nAtoms - 1, reductionLoops[threadIndex]);
        } else {
          execute(1, nAtoms, rescanTreeLoops[threadIndex]);
          execute(1, nAtoms, computeVolumeLoops[threadIndex]);
          execute(0, nAtoms - 1, reductionLoops[threadIndex]);
        }
      } catch (RuntimeException ex) {
        logger.warning("Runtime exception computing tree in thread: " + threadIndex);
        throw ex;
      } catch (Exception e) {
        String message = "Fatal exception computing tree in thread: " + threadIndex + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }

    /** Initialize the Overlap tree for a subset of the system. */
    private class ComputeTreeLoop extends IntegerForLoop {
      @Override
      public void run(int first, int last) throws Exception {
        // Compute the overlaps for a subset of atoms.
        int threadIndex = getThreadIndex();
        for (int slot = first; slot <= last; slot++) {
          localTree[threadIndex].computeAndAddChildrenR(slot);
        }
      }

      @Override
      public void start() {
        localTree[getThreadIndex()].initOverlapTree(
            coordinates, radii, volumes, gammas, ishydrogen);
      }
    }

    /** Compute the volume and derivatives for a subset of the system. */
    private class ComputeVolumeLoop extends IntegerForLoop {
      @Override
      public void run(int first, int last) throws Exception {
        int threadIndex = getThreadIndex();
        for (int slot = first; slot <= last; slot++) {
          double[] psi1i = new double[1];
          double[] f1i = new double[1];
          double[] p1i = new double[3];
          double[] psip1i = new double[1];
          double[] fp1i = new double[1];
          double[] pp1i = new double[3];
          double[] energy1i = new double[1];
          double[] fenergy1i = new double[1];
          double[] penergy1i = new double[3];
          localTree[threadIndex].computeVolumeUnderSlot2R(
              slot,
              psi1i,
              f1i,
              p1i,
              psip1i,
              fp1i,
              pp1i,
              energy1i,
              fenergy1i,
              penergy1i,
              threadIndex,
              grad,
              gradV,
              freeVolume,
              selfVolume);
          totalVolume.addAndGet(psi1i[0]);
          energy.addAndGet(energy1i[0]);
        }
      }
    }

    /** Rescan the tree based on updated radii and volumes for a subset of the system. */
    private class RescanTreeLoop extends IntegerForLoop {
      @Override
      public void run(int first, int last) throws Exception {
        // Compute the overlaps for a subset of atoms.
        int threadIndex = getThreadIndex();
        for (int slot = first; slot <= last; slot++) {
          localTree[threadIndex].rescanR(slot);
        }
      }

      @Override
      public void start() {
        localTree[getThreadIndex()].initRescanTreeV(
            coordinates, radii, volumes, gammas, ishydrogen);
      }
    }

    /** Reduce GaussVol gradient. */
    private class ReductionLoop extends IntegerForLoop {
      int threadID;

      @Override
      public void run(int lb, int ub) {
        grad.reduce(lb, ub);
        gradV.reduce(lb, ub);
        freeVolume.reduce(lb, ub);
        selfVolume.reduce(lb, ub);
        for (int i = lb; i <= ub; i++) {
          if (volumes[i] > 0) {
            gradV.set(0, i, gradV.get(i) / volumes[i]);
          }
        }
      }

      @Override
      public void start() {
        threadID = getThreadIndex();
      }
    }
  }
}
