// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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

import static ffx.potential.nonbonded.implicit.BornTanhRescaling.MAX_BORN_RADIUS;
import static ffx.potential.nonbonded.implicit.BornTanhRescaling.tanhRescaling;
import static ffx.potential.nonbonded.implicit.NeckIntegral.getNeckConstants;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.reduction.DoubleOp;
import edu.rit.pj.reduction.SharedDoubleArray;
import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;

/**
 * Parallel computation of Born radii via the Grycuk method.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class BornRadiiRegion extends ParallelRegion {

  private static final Logger logger = Logger.getLogger(BornRadiiRegion.class.getName());
  private static final double oneThird = 1.0 / 3.0;
  private static final double PI4_3 = 4.0 / 3.0 * PI;
  private static final double INVERSE_PI4_3 = 1.0 / PI4_3;
  private static final double PI_12 = PI / 12.0;
  private final BornRadiiLoop[] bornRadiiLoop;
  /** An ordered array of atoms in the system. */
  protected Atom[] atoms;
  /** Periodic boundary conditions and symmetry. */
  private Crystal crystal;
  /** Atomic coordinates for each symmetry operator. */
  private double[][][] sXYZ;
  /** Neighbor lists for each atom and symmetry operator. */
  private int[][][] neighborLists;
  /** Base radius of each atom. */
  private double[] baseRadius;
  /** Descreen radius of each atom. */
  private double[] descreenRadius;
  /**
   * Overlap scale factor for each atom, when using the Hawkins, Cramer & Truhlar pairwise
   * descreening algorithm.
   *
   * <p>G. D. Hawkins, C. J. Cramer and D. G. Truhlar, "Parametrized Models of Aqueous Free Energies
   * of Solvation Based on Pairwise Descreening of Solute Atomic Charges from a Dielectric Medium",
   * J. Phys. Chem., 100, 19824-19839 (1996).
   */
  private double[] overlapScale;
  /**
   * Sneck scaling parameter for each atom. Set based on maximum Sneck scaling parameter and number
   * of bound non-hydrogen atoms
   */
  private double[] neckScale;
  /**
   * Perfect born radius values
   */
  private final double[] perfectRadii;
  /**
   * Flag to turn on use of perfect Born radii.
   */
  private final boolean usePerfectRadii;
  /** Born radius of each atom. */
  private double[] born;
  /**
   * Flag to indicate if an atom should be included.
   */
  private boolean[] use;
  /** GK cut-off distance squared. */
  private double cut2;
  /**
   * Forces all atoms to be considered during Born radius updates.
   */
  private boolean nativeEnvironmentApproximation;
  /** If true, the descreening integral includes overlaps with the volume of the descreened atom */
  private final boolean perfectHCTScale;
  private double descreenOffset = 0.0;
  /**
   * The Born radius for each atom.
   */
  private SharedDoubleArray sharedBorn;
  /**
   * Boolean indicating whether or not to print all Born radii for an input molecular system. This is
   * turned off after the first round of printing.
   */
  private boolean verboseRadii;
  /**
   * Boolean indicating whether or not to use the neck volume correction for the implicit solvent
   */
  private final boolean neckCorrection;
  /**
   * Boolean indicating whether or not to use the tanh volume correction for the implicit solvent
   */
  private final boolean tanhCorrection;
  /**
   * This is the Born Integral prior to rescaling with a tanh function, which is a quantity needed
   * for the computing the derivative of the energy with respect to atomic coordinates.
   */
  private double[] unscaledBornIntegral;

  /**
   * BornRadiiRegion Constructor.
   *
   * @param nt Number of threads.
   * @param nAtoms Number of atoms.
   * @param forceField The ForceField in use.
   * @param neckCorrection Perform a neck correction.
   * @param tanhCorrection Perform a tanh correction.
   * @param perfectHCTScale Use "perfect" HCT scale factors.
   */
  public BornRadiiRegion(int nt, int nAtoms, ForceField forceField, boolean neckCorrection,
      boolean tanhCorrection, boolean perfectHCTScale) {
    bornRadiiLoop = new BornRadiiLoop[nt];
    for (int i = 0; i < nt; i++) {
      bornRadiiLoop[i] = new BornRadiiLoop();
    }
    verboseRadii = forceField.getBoolean("VERBOSE_BORN_RADII", false);
    this.perfectHCTScale = perfectHCTScale;
    this.neckCorrection = neckCorrection;
    this.tanhCorrection = tanhCorrection;
    if (verboseRadii && logger.isLoggable(Level.FINER)) {
      logger.finer(" Verbose Born radii.");
    }

    if (tanhCorrection) {
      unscaledBornIntegral = new double[nAtoms];
    }

    CompositeConfiguration compositeConfiguration = forceField.getProperties();
    usePerfectRadii = forceField.getBoolean("PERFECT_RADII", false);
    String[] radii = compositeConfiguration.getStringArray("perfect-radius");
    if (usePerfectRadii) {
      perfectRadii = new double[nAtoms];
      if (radii != null && radii.length > 0) {
        if (logger.isLoggable(Level.FINER)) {
          logger.finer(format(" Reading %d perfect-radius records .", radii.length));
        }
        for (String radius : radii) {
          String[] tokens = radius.trim().split(" +");
          if (tokens.length == 2) {
            // Input records should be from 1 to the number of atoms (subtract 1 to index from 0).
            int index = Integer.parseInt(tokens[0]) - 1;
            double value = Double.parseDouble(tokens[1]);
            if (logger.isLoggable(Level.FINER)) {
              logger.finer(format(" perfect-radius %d %16.8f", index, value));
            }
            if (index >= 0 && index < nAtoms && value > 0.0) {
              perfectRadii[index] = value;
            }
          } else {
            logger.warning(format(" Could not parse perfect-radius line %s", radius));
          }
        }
      }
    } else {
      perfectRadii = null;
    }
  }

  /**
   * Return perfect Born radii read in as keywords, or base radii if perfect radii are not available.
   *
   * @return Array of perfect Born radii.
   */
  public double[] getPerfectRadii() {
    int nAtoms = atoms.length;

    // Start with base radii.
    double[] radii = new double[nAtoms];
    arraycopy(baseRadius, 0, radii, 0, nAtoms);

    // Load perfect radii and return.
    if (usePerfectRadii) {
      for (int i = 0; i < nAtoms; i++) {
        double perfectRadius = perfectRadii[i];
        if (perfectRadius > 0.0) {
          radii[i] = perfectRadius;
        }
      }
    }

    return radii;
  }

  @Override
  public void finish() {
    int nAtoms = atoms.length;

    // Load perfect radii and return.
    if (usePerfectRadii) {
      for (int i = 0; i < nAtoms; i++) {
        double perfectRadius = perfectRadii[i];
        if (!use[i] || perfectRadius <= 0.0) {
          born[i] = baseRadius[i];
        } else {
          born[i] = perfectRadius;
        }
      }
      return;
    }

    for (int i = 0; i < nAtoms; i++) {
      final double baseRi = baseRadius[i];
      if (!use[i]) {
        born[i] = baseRi;
      } else {
        // A positive integral of 1/r^6 over the solute outside atom i.
        double soluteIntegral = -sharedBorn.get(i);
        if (tanhCorrection) {
          // Scale up the integral to account for interstitial spaces.
          unscaledBornIntegral[i] = soluteIntegral;
          soluteIntegral = tanhRescaling(soluteIntegral, baseRi);
        }
        // The total integral assumes no solute outside atom i, then subtracts away solute descreening.
        double sum = PI4_3 / (baseRi * baseRi * baseRi) - soluteIntegral;
        // Due to solute atomic overlaps, in rare cases the sum can be less than zero.
        if (sum <= 0.0) {
          born[i] = MAX_BORN_RADIUS;
          if (verboseRadii) {
            logger.info(format(
                " Born Integral < 0 for atom %d; set Born radius to %12.6f (Base Radius: %12.6f)",
                i + 1, born[i], baseRadius[i]));
          }
        } else {
          born[i] = pow(INVERSE_PI4_3 * sum, -oneThird);
          if (born[i] < baseRi) {
            born[i] = baseRi;
            if (verboseRadii) {
              logger.info(
                  format(" Born radius < Base Radius for atom %d: set Born radius to %12.6f", i + 1,
                      baseRi));
            }
          } else if (born[i] > MAX_BORN_RADIUS) {
            born[i] = MAX_BORN_RADIUS;
            if (verboseRadii) {
              logger.info(
                  format(" Born radius > 50.0 Angstroms for atom %d: set Born radius to %12.6f",
                      i + 1, baseRi));
            }
          } else if (isInfinite(born[i]) || isNaN(born[i])) {
            born[i] = baseRi;
            if (verboseRadii) {
              logger.info(
                  format(" Born radius NaN / Infinite for atom %d; set Born radius to %12.6f", i + 1,
                      baseRi));
            }
          } else {
            if (verboseRadii) {
              logger.info(
                  format(" Set Born radius for atom %d to %12.6f (Base Radius: %2.6f)", i + 1,
                      born[i], baseRi));
            }
          }
        }
      }
    }

    if (verboseRadii) {
      // Only log the Born radii once.
      logger.info(" Disabling verbose radii printing.");
      verboseRadii = false;
    }

  }

  public void init(
      Atom[] atoms,
      Crystal crystal,
      double[][][] sXYZ,
      int[][][] neighborLists,
      double[] baseRadius,
      double[] descreenRadius,
      double[] overlapScale,
      double[] neckScale,
      double descreenOffset,
      boolean[] use,
      double cut2,
      boolean nativeEnvironmentApproximation,
      double[] born) {
    this.atoms = atoms;
    this.crystal = crystal;
    this.sXYZ = sXYZ;
    this.neighborLists = neighborLists;
    this.baseRadius = baseRadius;
    this.descreenRadius = descreenRadius;
    this.overlapScale = overlapScale;
    this.neckScale = neckScale;
    this.descreenOffset = descreenOffset;
    this.use = use;
    this.cut2 = cut2;
    this.nativeEnvironmentApproximation = nativeEnvironmentApproximation;
    this.born = born;
  }

  @Override
  public void run() {
    if (!usePerfectRadii) {
      try {
        int nAtoms = atoms.length;
        execute(0, nAtoms - 1, bornRadiiLoop[getThreadIndex()]);
      } catch (Exception e) {
        String message = "Fatal exception computing Born radii in thread " + getThreadIndex() + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
  }

  @Override
  public void start() {
    int nAtoms = atoms.length;
    if (sharedBorn == null || sharedBorn.length() < nAtoms) {
      sharedBorn = new SharedDoubleArray(nAtoms);
    }
    for (int i = 0; i < nAtoms; i++) {
      sharedBorn.set(i, 0.0);
    }
  }

  public double[] getBorn() {
    return born;
  }

  public double[] getUnscaledBornIntegral() {
    return unscaledBornIntegral;
  }

  /**
   * Compute Born radii for a range of atoms via the Grycuk method.
   *
   * @since 1.0
   */
  private class BornRadiiLoop extends IntegerForLoop {

    private double[] localBorn;

    @Override
    public void finish() {
      sharedBorn.reduce(localBorn, DoubleOp.SUM);
    }

    @Override
    public void run(int lb, int ub) {
      int nSymm = crystal.spaceGroup.symOps.size();
      if (nSymm == 0) {
        nSymm = 1;
      }
      double[] x = sXYZ[0][0];
      double[] y = sXYZ[0][1];
      double[] z = sXYZ[0][2];
      for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
        double[][] xyz = sXYZ[iSymOp];
        for (int i = lb; i <= ub; i++) {
          if (!nativeEnvironmentApproximation && !use[i]) {
            continue;
          }
          final double integralStartI = max(baseRadius[i], descreenRadius[i]) + descreenOffset;
          final double descreenRi = descreenRadius[i];
          final double xi = x[i];
          final double yi = y[i];
          final double zi = z[i];
          int[] list = neighborLists[iSymOp][i];
          for (int k : list) {
            final double integralStartK = max(baseRadius[k], descreenRadius[k]) + descreenOffset;
            final double descreenRk = descreenRadius[k];
            if (!nativeEnvironmentApproximation && !use[k]) {
              continue;
            }

            // No necks will be computed unless the overlapScale is greater than 0.0 (e.g., for hydrogen).
            double mixedNeckScale = 0.5 * (neckScale[i] + neckScale[k]);

            if (i != k) {
              final double xr = xyz[0][k] - xi;
              final double yr = xyz[1][k] - yi;
              final double zr = xyz[2][k] - zi;
              final double r2 = crystal.image(xr, yr, zr);
              if (r2 > cut2) {
                continue;
              }
              final double r = sqrt(r2);
              // Atom i being descreeened by atom k.
              double sk = overlapScale[k];
              // Non-descreening atoms (such as hydrogen) will have an sk of 0.0
              if (sk > 0.0) {
                double descreenIK = descreen(r, r2, integralStartI, descreenRk, sk);
                localBorn[i] += descreenIK;
                if (neckCorrection) {
                  localBorn[i] += neckDescreen(r, integralStartI, descreenRk, mixedNeckScale);
                }
              }

              // Atom k being descreeened by atom i.
              double si = overlapScale[i];
              if (si > 0.0) {
                double descreenKI = descreen(r, r2, integralStartK, descreenRi, si);
                localBorn[k] += descreenKI;
                if (neckCorrection) {
                  localBorn[k] += neckDescreen(r, integralStartK, descreenRi, mixedNeckScale);
                }
              }

            } else if (iSymOp > 0) {
              final double xr = xyz[0][k] - xi;
              final double yr = xyz[1][k] - yi;
              final double zr = xyz[2][k] - zi;
              final double r2 = crystal.image(xr, yr, zr);
              if (r2 > cut2) {
                continue;
              }
              final double r = sqrt(r2);
              // Atom i being descreeened by atom k.
              double sk = overlapScale[k];
              if (sk > 0.0) {
                localBorn[i] += descreen(r, r2, integralStartI, descreenRk, sk);
                if (neckCorrection) {
                  localBorn[i] += neckDescreen(r, integralStartI, descreenRk, mixedNeckScale);
                }
              }
              // For symmetry mates, atom k is not descreeened by atom i.
            }
          }
        }
      }
    }

    @Override
    public void start() {
      int nAtoms = atoms.length;
      if (localBorn == null || localBorn.length < nAtoms) {
        localBorn = new double[nAtoms];
      }
      fill(localBorn, 0.0);
    }

    /**
     * Compute the integral of 1/r^6 over a neck region.
     *
     * @param r atomic separation.
     * @param radius base radius of the atom being descreened.
     * @param radiusK radius of the atom doing the descreening.
     * @param sneck Sneck scaling factor, scaled based on number of bound non-hydrogen atoms.
     * @return this contribution to the descreening integral.
     */
    private double neckDescreen(double r, double radius, double radiusK, double sneck) {
      double radiusWater = 1.4;

      // If atoms are too widely separated there is no neck formed.
      if (r > radius + radiusK + 2.0 * radiusWater) {
        return 0.0;
      }

      // Get Aij and Bij based on parameterization by Corrigan et al.
      double[] constants = getNeckConstants(radius, radiusK);

      double Aij = constants[0];
      double Bij = constants[1];
      //logger.info(format("rhoi %2.4f rhoj %2.4f separation %2.4f Aij %2.10f Bij %2.2f",radius,radiusK,r,Aij,Bij));

      double rMinusBij = r - Bij;
      double radiiMinusr = radius + radiusK + 2.0 * radiusWater - r;
      double power1 = rMinusBij * rMinusBij * rMinusBij * rMinusBij;
      double power2 = radiiMinusr * radiiMinusr * radiiMinusr * radiiMinusr;

      // Use Aij and Bij to get neck integral using Equations 13 and 14 from Aguilar/Onufriev 2010 paper
      // Sneck may be based on the number of heavy atoms bound to the atom being descreened.
      double neckIntegral = PI4_3 * sneck * Aij * power1 * power2;

      return -neckIntegral;
    }

    private double descreen(double r, double r2, double radius, double radiusK, double hctScale) {
      if (perfectHCTScale) {
        return perfectHCTIntegral(r, r2, radius, radiusK, hctScale);
      } else {
        return integral(r, r2, radius, radiusK * hctScale);
      }
    }

    /**
     * Use pairwise descreening to compute integral of 1/r^6.
     *
     * @param r atomic separation.
     * @param r2 atomic separation squared.
     * @param radius base radius of the atom being descreened.
     * @param scaledRadius scaled radius of the atom doing the descreening.
     * @return this contribution to the descreening integral.
     */
    private double integral(double r, double r2, double radius, double scaledRadius) {
      double integral = 0.0;
      // Descreen only if the scaledRadius is greater than zero.
      // and atom I does not engulf atom K.
      if (scaledRadius > 0.0 && (radius < r + scaledRadius)) {
        // Atom i is engulfed by atom k.
        if (radius + r < scaledRadius) {
          final double upper = scaledRadius - r;
          integral = (PI4_3 * (1.0 / (upper * upper * upper) - 1.0 / (radius * radius * radius)));
        }

        // Upper integration bound is always the same.
        double upper = r + scaledRadius;

        // Lower integration bound depends on atoms sizes and separation.
        double lower;
        if (radius + r < scaledRadius) {
          // Atom i is engulfed by atom k.
          lower = scaledRadius - r;
        } else if (r < radius + scaledRadius) {
          // Atoms are overlapped, begin integration from ri.
          lower = radius;
        } else {
          // No overlap between atoms.
          lower = r - scaledRadius;
        }

        double l2 = lower * lower;
        double l4 = l2 * l2;
        double lr = lower * r;
        double l4r = l4 * r;
        double u2 = upper * upper;
        double u4 = u2 * u2;
        double ur = upper * r;
        double u4r = u4 * r;
        double scaledRk2 = scaledRadius * scaledRadius;
        double term =
            (3.0 * (r2 - scaledRk2) + 6.0 * u2 - 8.0 * ur) / u4r
                - (3.0 * (r2 - scaledRk2) + 6.0 * l2 - 8.0 * lr) / l4r;
        integral -= PI_12 * term;
      }
      return integral;
    }

    private double perfectHCTIntegral(double r, double r2,
        double radius, double radiusK, double perfectHCT) {
      double integral = 0.0;
      // Descreen only if the scaledRadius is greater than zero.
      // and atom I does not engulf atom K.
      if (radiusK > 0.0 && (radius < r + radiusK)) {
        // Atom i is engulfed by atom k.
        // TODO: fix double counting of overlaps
        if (radius + r < radiusK) {
          final double upper = radiusK - r;
          integral = (PI4_3 * (1.0 / (upper * upper * upper) - 1.0 / (radius * radius * radius)));
        }

        // Upper integration bound is always the same.
        double upper = r + radiusK;

        // Lower integration bound depends on atoms sizes and separation.
        double lower;
        if (radius + r < radiusK) {
          // Atom i is engulfed by atom k.
          lower = radiusK - r;
        } else if (r < radius + radiusK) {
          // Atoms are overlapped, begin integration from ri.
          // TODO: fix double counting of the overlap
          lower = radius;
        } else {
          // No overlap between atoms.
          lower = r - radiusK;
        }

        double l2 = lower * lower;
        double l4 = l2 * l2;
        double lr = lower * r;
        double l4r = l4 * r;
        double u2 = upper * upper;
        double u4 = u2 * u2;
        double ur = upper * r;
        double u4r = u4 * r;
        double scaledK2 = radiusK * radiusK;
        double term =
            (3.0 * (r2 - scaledK2) + 6.0 * u2 - 8.0 * ur) / u4r
                - (3.0 * (r2 - scaledK2) + 6.0 * l2 - 8.0 * lr) / l4r;
        integral -= PI_12 * term;
      }

      return perfectHCT * integral;
    }
  }
}
