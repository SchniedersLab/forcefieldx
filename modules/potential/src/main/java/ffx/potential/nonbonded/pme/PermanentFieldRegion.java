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
package ffx.potential.nonbonded.pme;

import static ffx.numerics.special.Erf.erfc;
import static ffx.potential.parameters.MultipoleType.t000;
import static ffx.potential.parameters.MultipoleType.t001;
import static ffx.potential.parameters.MultipoleType.t002;
import static ffx.potential.parameters.MultipoleType.t010;
import static ffx.potential.parameters.MultipoleType.t011;
import static ffx.potential.parameters.MultipoleType.t020;
import static ffx.potential.parameters.MultipoleType.t100;
import static ffx.potential.parameters.MultipoleType.t101;
import static ffx.potential.parameters.MultipoleType.t110;
import static ffx.potential.parameters.MultipoleType.t200;
import static ffx.utilities.Constants.ELEC_ANG_TO_DEBYE;
import static java.lang.String.format;
import static java.util.Arrays.copyOf;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedInteger;
import edu.rit.util.Range;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.MaskingInterface;
import ffx.potential.nonbonded.NeighborList;
import ffx.potential.nonbonded.ReciprocalSpace;
import ffx.potential.parameters.ForceField;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Parallel computation of the permanent field.
 *
 * <p>This class can be executed by a ParallelTeam with exactly 2 threads.
 *
 * <p>The Real Space and Reciprocal Space Sections will be run concurrently, each with the number of
 * threads defined by their respective ParallelTeam instances.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PermanentFieldRegion extends ParallelRegion implements MaskingInterface {

  private static final Logger logger = Logger.getLogger(PermanentFieldRegion.class.getName());
  /**
   * Constant applied to multipole interactions.
   */
  private static final double oneThird = 1.0 / 3.0;
  /**
   * Specify inter-molecular softcore.
   */
  private final boolean intermolecularSoftcore;
  /**
   * Specify intra-molecular softcore.
   */
  private final boolean intramolecularSoftcore;
  /**
   * Dimensions of [nsymm][nAtoms][3]
   */
  public double[][][] inducedDipole;
  public double[][][] inducedDipoleCR;
  /**
   * Polarization groups.
   */
  protected int[][] ip11;
  /**
   * An ordered array of atoms in the system.
   */
  private Atom[] atoms;
  /**
   * Unit cell and spacegroup information.
   */
  private Crystal crystal;
  /**
   * Dimensions of [nsymm][xyz][nAtoms].
   */
  private double[][][] coordinates;
  /**
   * Dimensions of [nsymm][nAtoms][10]
   */
  private double[][][] globalMultipole;
  /**
   * Neighbor lists, including atoms beyond the real space cutoff. [nsymm][nAtoms][nAllNeighbors]
   */
  private int[][][] neighborLists;
  /**
   * Neighbor list cells for Octree method
   */
  private NeighborList.Cell[][][] cells;
  /**
   * Neighbor lists, without atoms beyond the preconditioner cutoff.
   * [nSymm][nAtoms][nIncludedNeighbors]
   */
  private int[][][] preconditionerLists;
  /**
   * Number of neighboring atoms within the preconditioner cutoff. [nSymm][nAtoms]
   */
  private int[][] preconditionerCounts;
  /**
   * When computing the polarization energy at Lambda there are 3 pieces.
   *
   * <p>1.) Upol(1) = The polarization energy computed normally (ie. system with ligand).
   *
   * <p>2.) Uenv = The polarization energy of the system without the ligand.
   *
   * <p>3.) Uligand = The polarization energy of the ligand by itself.
   *
   * <p>Upol(L) = L*Upol(1) + (1-L)*(Uenv + Uligand)
   *
   * <p>Set the "use" array to true for all atoms for part 1. Set the "use" array to true for all
   * atoms except the ligand for part 2. Set the "use" array to true only for the ligand atoms for
   * part 3.
   *
   * <p>The "use" array can also be employed to turn off atoms for computing the electrostatic
   * energy of sub-structures.
   */
  private boolean[] use;
  /**
   * Molecule number for each atom.
   */
  private int[] molecule;
  private double[] ipdamp;
  private double[] thole;
  /**
   * Masking of 1-2, 1-3 and 1-4 interactions.
   */
  private int[][] mask12;
  private int[][] mask13;
  private int[][] mask14;
  /**
   * The current LambdaMode of this PME instance (or OFF for no lambda dependence).
   */
  private LambdaMode lambdaMode = LambdaMode.OFF;
  /**
   * Reciprocal space instance.
   */
  private ReciprocalSpace reciprocalSpace;
  private boolean reciprocalSpaceTerm;
  private double off2;
  private double preconditionerCutoff;
  private double an0, an1, an2;
  private double aewald;
  /**
   * Neighbor lists, without atoms beyond the real space cutoff. [nSymm][nAtoms][nIncludedNeighbors]
   */
  private int[][][] realSpaceLists;
  /**
   * Number of neighboring atoms within the real space cutoff. [nSymm][nAtoms]
   */
  private int[][] realSpaceCounts;
  private Range[] realSpaceRanges;
  private IntegerSchedule permanentSchedule;
  /**
   * Field array.
   */
  private AtomicDoubleArray3D field;
  /**
   * Chain rule field array.
   */
  private AtomicDoubleArray3D fieldCR;
  private ScaleParameters scaleParameters;
  private final PermanentRealSpaceFieldSection permanentRealSpaceFieldSection;
  private final PermanentReciprocalSection permanentReciprocalSection;
  /**
   * Timing variables.
   */
  private PMETimings pmeTimings;

  public PermanentFieldRegion(ParallelTeam pt, ForceField forceField, boolean lambdaTerm) {
    permanentRealSpaceFieldSection = new PermanentRealSpaceFieldSection(pt);
    permanentReciprocalSection = new PermanentReciprocalSection();

    // Flag to indicate application of an intermolecular softcore potential.
    if (lambdaTerm) {
      intermolecularSoftcore = forceField.getBoolean("INTERMOLECULAR_SOFTCORE", false);
      intramolecularSoftcore = forceField.getBoolean("INTRAMOLECULAR_SOFTCORE", false);
    } else {
      intermolecularSoftcore = false;
      intramolecularSoftcore = false;
    }
  }

  /**
   * Apply permanent field masking rules.
   *
   * @param i     The atom whose masking rules should be applied.
   * @param is14  True if atom i and the current atom are 1-4 to each other.
   * @param masks One or more masking arrays.
   */
  @Override
  public void applyMask(int i, boolean[] is14, double[]... masks) {
    if (ip11[i] != null) {
      double[] energyMask = masks[0];
      var m12 = mask12[i];
      for (int value : m12) {
        energyMask[value] = scaleParameters.p12scale;
      }
      var m13 = mask13[i];
      for (int value : m13) {
        energyMask[value] = scaleParameters.p13scale;
      }
      var m14 = mask14[i];
      for (int value : m14) {
        energyMask[value] = scaleParameters.p14scale;
        for (int k : ip11[i]) {
          if (k == value) {
            energyMask[value] = scaleParameters.intra14Scale * scaleParameters.p14scale;
            break;
          }
        }
      }
      // Apply group based polarization masking rule.
      double[] inductionMask = masks[1];
      for (int index : ip11[i]) {
        inductionMask[index] = scaleParameters.d11scale;
      }
    }
  }

  public void initTimings() {
    permanentRealSpaceFieldSection.initTimings();
  }

  public long getRealSpacePermTime() {
    return permanentRealSpaceFieldSection.time;
  }

  public long getInitTime(int threadId) {
    return permanentRealSpaceFieldSection.permanentRealSpaceFieldRegion.initializationLoop[threadId].time;
  }

  public long getPermTime(int threadId) {
    return permanentRealSpaceFieldSection.permanentRealSpaceFieldRegion.permanentRealSpaceFieldLoop[threadId].time;
  }

  public void init(
      Atom[] atoms,
      Crystal crystal,
      double[][][] coordinates,
      double[][][] globalMultipole,
      double[][][] inducedDipole,
      double[][][] inducedDipoleCR,
      int[][][] neighborLists,
      NeighborList.Cell[][][] cells,
      ScaleParameters scaleParameters,
      boolean[] use,
      int[] molecule,
      double[] ipdamp,
      double[] thole,
      int[][] ip11,
      int[][] mask12,
      int[][] mask13,
      int[][] mask14,
      LambdaMode lambdaMode,
      boolean reciprocalSpaceTerm,
      ReciprocalSpace reciprocalSpace,
      EwaldParameters ewaldParameters,
      PCGSolver pcgSolver,
      IntegerSchedule permanentSchedule,
      RealSpaceNeighborParameters realSpaceNeighborParameters,
      AtomicDoubleArray3D field,
      AtomicDoubleArray3D fieldCR) {
    this.atoms = atoms;
    this.crystal = crystal;
    this.coordinates = coordinates;
    this.globalMultipole = globalMultipole;
    this.inducedDipole = inducedDipole;
    this.inducedDipoleCR = inducedDipoleCR;
    this.neighborLists = neighborLists;
    this.cells = cells;
    this.scaleParameters = scaleParameters;
    this.use = use;
    this.molecule = molecule;
    this.ipdamp = ipdamp;
    this.thole = thole;
    this.ip11 = ip11;
    this.mask12 = mask12;
    this.mask13 = mask13;
    this.mask14 = mask14;
    this.lambdaMode = lambdaMode;
    this.reciprocalSpaceTerm = reciprocalSpaceTerm;
    this.reciprocalSpace = reciprocalSpace;
    if (pcgSolver != null) {
      this.preconditionerCutoff = pcgSolver.getPreconditionerCutoff();
      this.preconditionerLists = pcgSolver.getPreconditionerLists();
      this.preconditionerCounts = pcgSolver.getPreconditionerCounts();
    }
    this.aewald = ewaldParameters.aewald;
    this.an0 = ewaldParameters.an0;
    this.an1 = ewaldParameters.an1;
    this.an2 = ewaldParameters.an2;
    this.off2 = ewaldParameters.off2;
    this.permanentSchedule = permanentSchedule;
    this.realSpaceLists = realSpaceNeighborParameters.realSpaceLists;
    this.realSpaceCounts = realSpaceNeighborParameters.realSpaceCounts;
    this.realSpaceRanges = realSpaceNeighborParameters.realSpaceRanges;
    this.field = field;
    this.fieldCR = fieldCR;
  }

  /**
   * Remove permanent field masking rules.
   *
   * @param i     The atom whose masking rules should be removed.
   * @param is14  True if atom i and the current atom are 1-4 to each other.
   * @param masks One or more masking arrays.
   */
  @Override
  public void removeMask(int i, boolean[] is14, double[]... masks) {
    if (ip11[i] != null) {
      double[] energyMask = masks[0];
      var m12 = mask12[i];
      for (int value : m12) {
        energyMask[value] = 1.0;
      }
      var m13 = mask13[i];
      for (int value : m13) {
        energyMask[value] = 1.0;
      }
      var m14 = mask14[i];
      for (int value : m14) {
        energyMask[value] = 1.0;
        for (int k : ip11[i]) {
          if (k == value) {
            energyMask[value] = 1.0;
            break;
          }
        }
      }
      // Apply group based polarization masking rule.
      double[] inductionMask = masks[1];
      for (int index : ip11[i]) {
        inductionMask[index] = 1.0;
      }
    }
  }

  @Override
  public void run() {
    try {
      execute(permanentRealSpaceFieldSection, permanentReciprocalSection);
    } catch (RuntimeException e) {
      String message = "Runtime exception computing the permanent multipole field.\n";
      logger.log(Level.WARNING, message, e);
      throw e;
    } catch (Exception e) {
      String message = "Fatal exception computing the permanent multipole field.\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  /**
   * Computes the Permanent Multipole Real Space Field.
   */
  private class PermanentRealSpaceFieldSection extends ParallelSection {

    private final PermanentRealSpaceFieldRegion permanentRealSpaceFieldRegion;
    private final ParallelTeam parallelTeam;
    protected long time;

    PermanentRealSpaceFieldSection(ParallelTeam pt) {
      this.parallelTeam = pt;
      int nt = pt.getThreadCount();
      permanentRealSpaceFieldRegion = new PermanentRealSpaceFieldRegion(nt);
    }

    public void initTimings() {
      time = 0;
      permanentRealSpaceFieldRegion.initTimings();
    }

    @Override
    public void run() {
      try {
        time -= System.nanoTime();
        parallelTeam.execute(permanentRealSpaceFieldRegion);
        time += System.nanoTime();
      } catch (RuntimeException e) {
        String message = "Fatal exception computing the real space field.\n";
        logger.log(Level.WARNING, message, e);
      } catch (Exception e) {
        String message = "Fatal exception computing the real space field.\n";
        logger.log(Level.SEVERE, message, e);
      }
    }
  }

  /**
   * Compute the permanent multipole reciprocal space contribution to the electric potential, field,
   * etc. using the number of threads specified by the ParallelTeam used to construct the
   * ReciprocalSpace instance.
   */
  private class PermanentReciprocalSection extends ParallelSection {

    @Override
    public void run() {
      if (reciprocalSpaceTerm && aewald > 0.0) {
        reciprocalSpace.performConvolution();
      }
    }
  }

  private class PermanentRealSpaceFieldRegion extends ParallelRegion {

    private final InitializationLoop[] initializationLoop;
    private final PermanentRealSpaceFieldLoop[] permanentRealSpaceFieldLoop;
    private final SharedInteger sharedCount;
    private final int threadCount;

    PermanentRealSpaceFieldRegion(int nt) {
      threadCount = nt;
      initializationLoop = new InitializationLoop[threadCount];
      permanentRealSpaceFieldLoop = new PermanentRealSpaceFieldLoop[threadCount];
      sharedCount = new SharedInteger();
      for (int i = 0; i < threadCount; i++) {
        initializationLoop[i] = new InitializationLoop();
        permanentRealSpaceFieldLoop[i] = new PermanentRealSpaceFieldLoop();
      }
    }

    public void initTimings() {
      for (int i = 0; i < threadCount; i++) {
        permanentRealSpaceFieldLoop[i].time = 0;
        initializationLoop[i].time = 0;
      }
    }

    @Override
    public void finish() {
      if (realSpaceRanges == null) {
        logger.severe(" RealSpaceRange array is null");
      }

      int nAtoms = atoms.length;

      // Load balancing.
      int id = 0;
      int goal = sharedCount.get() / threadCount;
      int num = 0;
      int start = 0;

      for (int i = 0; i < nAtoms; i++) {
        List<SymOp> symOps = crystal.spaceGroup.symOps;
        int nSymm = symOps.size();
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
          num += realSpaceCounts[iSymm][i];
        }
        if (num >= goal) {
          // Last thread gets the remaining atoms in its range.
          if (id == threadCount - 1) {
            realSpaceRanges[id] = new Range(start, nAtoms - 1);
            break;
          }

          realSpaceRanges[id] = new Range(start, i);

          // Reset the count.
          num = 0;

          // Next thread.
          id++;

          // Next range starts at i+1.
          start = i + 1;

          // Out of atoms. Threads remaining get a null range.
          if (start == nAtoms) {
            for (int j = id; j < threadCount; j++) {
              realSpaceRanges[j] = null;
            }
            break;
          }
        } else if (i == nAtoms - 1) {

          // Last atom without reaching goal for current thread.
          realSpaceRanges[id] = new Range(start, nAtoms - 1);
          for (int j = id + 1; j < threadCount; j++) {
            realSpaceRanges[j] = null;
          }
        }
      }
    }

    @Override
    public void run() {
      int threadIndex = getThreadIndex();
      try {
        int nAtoms = atoms.length;
        execute(0, nAtoms - 1, initializationLoop[threadIndex]);
        execute(0, nAtoms - 1, permanentRealSpaceFieldLoop[threadIndex]);
      } catch (RuntimeException e) {
        String message = "Runtime exception computing the real space field.\n";
        logger.log(Level.SEVERE, message, e);
      } catch (Exception e) {
        String message =
            "Fatal exception computing the real space field in thread " + getThreadIndex() + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }

    @Override
    public void start() {
      sharedCount.set(0);
    }

    private class InitializationLoop extends IntegerForLoop {

      protected long time;

      @Override
      public void start() {
        time -= System.nanoTime();
      }

      @Override
      public void finish() {
        time += System.nanoTime();
      }

      @Override
      public void run(int lb, int ub) {
        // Initialize the induced dipole arrays.
        List<SymOp> symOps = crystal.spaceGroup.symOps;
        int nSymm = symOps.size();
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
          double[][] ind0 = inducedDipole[0];
          double[][] indCR0 = inducedDipoleCR[0];
          for (int i = lb; i <= ub; i++) {
            double[] ind = ind0[i];
            double[] indCR = indCR0[i];
            ind[0] = 0.0;
            ind[1] = 0.0;
            ind[2] = 0.0;
            indCR[0] = 0.0;
            indCR[1] = 0.0;
            indCR[2] = 0.0;
          }
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return IntegerSchedule.fixed();
      }

    }

    private class PermanentRealSpaceFieldLoop extends IntegerForLoop {

      protected long time;
      private final double[] dx_local;
      private final double[][] transOp;
      private int threadID;
      private double[] inductionMaskLocal;
      private double[] energyMaskLocal;
      private int count;

      PermanentRealSpaceFieldLoop() {
        super();
        dx_local = new double[3];
        transOp = new double[3][3];
      }

      @Override
      public void finish() {
        sharedCount.addAndGet(count);
        time += System.nanoTime();
      }

      @Override
      public void run(int lb, int ub) {
        int[][] lists = neighborLists[0];
        int[][] ewalds = realSpaceLists[0];
        int[] counts = realSpaceCounts[0];
        int[][] preLists = preconditionerLists[0];
        int[] preCounts = preconditionerCounts[0];
        final double[] x = coordinates[0][0];
        final double[] y = coordinates[0][1];
        final double[] z = coordinates[0][2];
        final double[][] mpole = globalMultipole[0];
        // Loop over atom chunk.
        for (int i = lb; i <= ub; i++) {
          if (!use[i]) {
            continue;
          }

          // Zero out accumulation variables for the field at atom i.
          double fix = 0.0;
          double fiy = 0.0;
          double fiz = 0.0;
          double fixCR = 0.0;
          double fiyCR = 0.0;
          double fizCR = 0.0;

          final int moleculei = molecule[i];
          final double pdi = ipdamp[i];
          final double pti = thole[i];
          final double xi = x[i];
          final double yi = y[i];
          final double zi = z[i];
          final double[] globalMultipolei = mpole[i];
          final double ci = globalMultipolei[0];
          final double dix = globalMultipolei[t100];
          final double diy = globalMultipolei[t010];
          final double diz = globalMultipolei[t001];
          final double qixx = globalMultipolei[t200] * oneThird;
          final double qiyy = globalMultipolei[t020] * oneThird;
          final double qizz = globalMultipolei[t002] * oneThird;
          final double qixy = globalMultipolei[t110] * oneThird;
          final double qixz = globalMultipolei[t101] * oneThird;
          final double qiyz = globalMultipolei[t011] * oneThird;

          // Apply field masking rules.
          applyMask(i, null, energyMaskLocal, inductionMaskLocal);

          // Loop over the neighbor list.
          final int[] list = lists[i];
          counts[i] = 0;
          preCounts[i] = 0;
          int[] ewald = ewalds[i];
          int[] preList = preLists[i];
          for (int k : list) {
            if (!use[k]) {
              continue;
            }
            if (lambdaMode == LambdaMode.VAPOR) {
              boolean sameMolecule = (moleculei == molecule[k]);
              if ((intermolecularSoftcore && !sameMolecule)
                  || (intramolecularSoftcore && sameMolecule)) {
                continue;
              }
            }
            final double xk = x[k];
            final double yk = y[k];
            final double zk = z[k];
            dx_local[0] = xk - xi;
            dx_local[1] = yk - yi;
            dx_local[2] = zk - zi;
            final double r2 = crystal.image(dx_local);
            if (r2 <= off2) {
              count++;
              // Store a short neighbor list for the SCF.
              if (ewald.length <= counts[i]) {
                int len = ewald.length;
                ewalds[i] = copyOf(ewald, len + 10);
                ewald = ewalds[i];
              }
              ewald[counts[i]++] = k;
              final double xr = dx_local[0];
              final double yr = dx_local[1];
              final double zr = dx_local[2];
              final double pdk = ipdamp[k];
              final double ptk = thole[k];
              final double[] globalMultipolek = mpole[k];
              final double ck = globalMultipolek[t000];
              final double dkx = globalMultipolek[t100];
              final double dky = globalMultipolek[t010];
              final double dkz = globalMultipolek[t001];
              final double qkxx = globalMultipolek[t200] * oneThird;
              final double qkyy = globalMultipolek[t020] * oneThird;
              final double qkzz = globalMultipolek[t002] * oneThird;
              final double qkxy = globalMultipolek[t110] * oneThird;
              final double qkxz = globalMultipolek[t101] * oneThird;
              final double qkyz = globalMultipolek[t011] * oneThird;
              double r = sqrt(r2);

              // Store a short neighbor list for the SCF pre-conditioner.
              if (r < preconditionerCutoff) {
                if (preList.length <= preCounts[i]) {
                  int len = preList.length;
                  preLists[i] = copyOf(preList, len + 10);
                  preList = preLists[i];
                }
                preList[preCounts[i]++] = k;
              }

              // Calculate the error function damping terms.
              final double ralpha = aewald * r;
              final double exp2a = exp(-ralpha * ralpha);
              final double rr1 = 1.0 / r;
              final double rr2 = rr1 * rr1;
              final double bn0 = erfc(ralpha) * rr1;
              final double bn1 = (bn0 + an0 * exp2a) * rr2;
              final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
              final double bn3 = (5.0 * bn2 + an2 * exp2a) * rr2;

              // Compute the error function scaled and unscaled terms.
              double scale3 = 1.0;
              double scale5 = 1.0;
              double scale7 = 1.0;
              double damp = pdi * pdk;
              final double pgamma = min(pti, ptk);
              final double rdamp = r * damp;
              damp = -pgamma * rdamp * rdamp * rdamp;
              if (damp > -50.0) {
                double expdamp = exp(damp);
                scale3 = 1.0 - expdamp;
                scale5 = 1.0 - expdamp * (1.0 - damp);
                scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
              }
              final double scale = inductionMaskLocal[k];
              final double scalep = energyMaskLocal[k];
              // Thole damping multiplied by the group-based mask.
              final double dsc3 = scale3 * scale;
              final double dsc5 = scale5 * scale;
              final double dsc7 = scale7 * scale;
              // Thole damping multiplied by the energy mask.
              final double psc3 = scale3 * scalep;
              final double psc5 = scale5 * scalep;
              final double psc7 = scale7 * scalep;
              final double rr3 = rr1 * rr2;
              final double rr5 = 3.0 * rr3 * rr2;
              final double rr7 = 5.0 * rr5 * rr2;
              // 1.0 minus induction masks and Thole damping.
              final double drr3 = (1.0 - dsc3) * rr3;
              final double drr5 = (1.0 - dsc5) * rr5;
              final double drr7 = (1.0 - dsc7) * rr7;
              // 1.0 minus energy masks and Thole damping.
              final double prr3 = (1.0 - psc3) * rr3;
              final double prr5 = (1.0 - psc5) * rr5;
              final double prr7 = (1.0 - psc7) * rr7;
              final double dir = dix * xr + diy * yr + diz * zr;
              final double qix = 2.0 * (qixx * xr + qixy * yr + qixz * zr);
              final double qiy = 2.0 * (qixy * xr + qiyy * yr + qiyz * zr);
              final double qiz = 2.0 * (qixz * xr + qiyz * yr + qizz * zr);
              final double qir = (qix * xr + qiy * yr + qiz * zr) * 0.5;
              // Ewald field for atom k (no masking).
              final double bn123i = bn1 * ci + bn2 * dir + bn3 * qir;
              final double fkmx = xr * bn123i - bn1 * dix - bn2 * qix;
              final double fkmy = yr * bn123i - bn1 * diy - bn2 * qiy;
              final double fkmz = zr * bn123i - bn1 * diz - bn2 * qiz;
              // Correct Ewald field for over-counted induction interactions.
              final double ddr357i = drr3 * ci + drr5 * dir + drr7 * qir;
              final double fkdx = xr * ddr357i - drr3 * dix - drr5 * qix;
              final double fkdy = yr * ddr357i - drr3 * diy - drr5 * qiy;
              final double fkdz = zr * ddr357i - drr3 * diz - drr5 * qiz;
              field.add(threadID, k, fkmx - fkdx, fkmy - fkdy, fkmz - fkdz);
              // Correct Ewald field for over-counted energy interactions.
              final double prr357i = prr3 * ci + prr5 * dir + prr7 * qir;
              final double fkpx = xr * prr357i - prr3 * dix - prr5 * qix;
              final double fkpy = yr * prr357i - prr3 * diy - prr5 * qiy;
              final double fkpz = zr * prr357i - prr3 * diz - prr5 * qiz;
              fieldCR.add(threadID, k, fkmx - fkpx, fkmy - fkpy, fkmz - fkpz);
              final double dkr = dkx * xr + dky * yr + dkz * zr;
              final double qkx = 2.0 * (qkxx * xr + qkxy * yr + qkxz * zr);
              final double qky = 2.0 * (qkxy * xr + qkyy * yr + qkyz * zr);
              final double qkz = 2.0 * (qkxz * xr + qkyz * yr + qkzz * zr);
              final double qkr = (qkx * xr + qky * yr + qkz * zr) * 0.5;
              final double bn123k = bn1 * ck - bn2 * dkr + bn3 * qkr;
              // Ewald field for atom i (no masking).
              final double fimx = -xr * bn123k - bn1 * dkx + bn2 * qkx;
              final double fimy = -yr * bn123k - bn1 * dky + bn2 * qky;
              final double fimz = -zr * bn123k - bn1 * dkz + bn2 * qkz;
              // Correct Ewald field for over-counted induction interactions.
              final double drr357k = drr3 * ck - drr5 * dkr + drr7 * qkr;
              final double fidx = -xr * drr357k - drr3 * dkx + drr5 * qkx;
              final double fidy = -yr * drr357k - drr3 * dky + drr5 * qky;
              final double fidz = -zr * drr357k - drr3 * dkz + drr5 * qkz;
              fix += fimx - fidx;
              fiy += fimy - fidy;
              fiz += fimz - fidz;
              // Correct Ewald field for over-counted energy interactions.
              final double prr357k = prr3 * ck - prr5 * dkr + prr7 * qkr;
              final double fipx = -xr * prr357k - prr3 * dkx + prr5 * qkx;
              final double fipy = -yr * prr357k - prr3 * dky + prr5 * qky;
              final double fipz = -zr * prr357k - prr3 * dkz + prr5 * qkz;
              fixCR += fimx - fipx;
              fiyCR += fimy - fipy;
              fizCR += fimz - fipz;
            }
          }
          // Add in field contributions at Atom i.
          field.add(threadID, i, fix, fiy, fiz);
          fieldCR.add(threadID, i, fixCR, fiyCR, fizCR);
          // Remove field masking rules.
          removeMask(i, null, energyMaskLocal, inductionMaskLocal);
        }

        // Loop over symmetry mates.
        List<SymOp> symOps = crystal.spaceGroup.symOps;
        int nSymm = symOps.size();
        for (int iSymm = 1; iSymm < nSymm; iSymm++) {
          SymOp symOp = crystal.spaceGroup.getSymOp(iSymm);
          crystal.getTransformationOperator(symOp, transOp);
          lists = neighborLists[iSymm];
          ewalds = realSpaceLists[iSymm];
          counts = realSpaceCounts[iSymm];
          preLists = preconditionerLists[iSymm];
          preCounts = preconditionerCounts[iSymm];
          double[] xs = coordinates[iSymm][0];
          double[] ys = coordinates[iSymm][1];
          double[] zs = coordinates[iSymm][2];
          double[][] mpoles = globalMultipole[iSymm];

          // Loop over atoms in a chunk of the asymmetric unit.
          for (int i = lb; i <= ub; i++) {
            if (!use[i]) {
              continue;
            }
            // Zero out accumulation variables for the field at atom i.
            double fix = 0.0;
            double fiy = 0.0;
            double fiz = 0.0;
            final double pdi = ipdamp[i];
            final double pti = thole[i];
            final double[] multipolei = mpole[i];
            final double ci = multipolei[t000];
            final double dix = multipolei[t100];
            final double diy = multipolei[t010];
            final double diz = multipolei[t001];
            final double qixx = multipolei[t200] * oneThird;
            final double qiyy = multipolei[t020] * oneThird;
            final double qizz = multipolei[t002] * oneThird;
            final double qixy = multipolei[t110] * oneThird;
            final double qixz = multipolei[t101] * oneThird;
            final double qiyz = multipolei[t011] * oneThird;
            final double xi = x[i];
            final double yi = y[i];
            final double zi = z[i];

            // Loop over the neighbor list.
            final int[] list = lists[i];
            counts[i] = 0;
            preCounts[i] = 0;
            int[] ewald = ewalds[i];
            int[] preList = preLists[i];
            for (int k : list) {
              if (!use[k]) {
                continue;
              }
              final double xk = xs[k];
              final double yk = ys[k];
              final double zk = zs[k];
              dx_local[0] = xk - xi;
              dx_local[1] = yk - yi;
              dx_local[2] = zk - zi;
              final double r2 = crystal.image(dx_local);
              if (r2 <= off2) {
                count++;
                // Store a short neighbor list for the SCF.
                if (ewald.length <= counts[i]) {
                  int len = ewald.length;
                  ewalds[i] = copyOf(ewald, len + 10);
                  ewald = ewalds[i];
                }
                ewald[counts[i]++] = k;
                double selfScale = 1.0;
                if (i == k) {
                  selfScale = 0.5;
                }
                final double xr = dx_local[0];
                final double yr = dx_local[1];
                final double zr = dx_local[2];
                final double pdk = ipdamp[k];
                final double ptk = thole[k];
                final double[] multipolek = mpoles[k];
                final double ck = multipolek[t000];
                final double dkx = multipolek[t100];
                final double dky = multipolek[t010];
                final double dkz = multipolek[t001];
                final double qkxx = multipolek[t200] * oneThird;
                final double qkyy = multipolek[t020] * oneThird;
                final double qkzz = multipolek[t002] * oneThird;
                final double qkxy = multipolek[t110] * oneThird;
                final double qkxz = multipolek[t101] * oneThird;
                final double qkyz = multipolek[t011] * oneThird;
                final double r = sqrt(r2);
                if (r < preconditionerCutoff) {
                  if (preList.length <= preCounts[i]) {
                    int len = preList.length;
                    preLists[i] = copyOf(preList, len + 10);
                    preList = preLists[i];
                  }
                  preList[preCounts[i]++] = k;
                }

                // Calculate the error function damping terms.
                final double ralpha = aewald * r;
                final double exp2a = exp(-ralpha * ralpha);
                final double rr1 = 1.0 / r;
                final double rr2 = rr1 * rr1;
                final double bn0 = erfc(ralpha) * rr1;
                final double bn1 = (bn0 + an0 * exp2a) * rr2;
                final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                final double bn3 = (5.0 * bn2 + an2 * exp2a) * rr2;

                // Compute the error function scaled and unscaled terms.
                double scale3 = 1.0;
                double scale5 = 1.0;
                double scale7 = 1.0;
                double damp = pdi * pdk;
                final double pgamma = min(pti, ptk);
                final double rdamp = r * damp;
                damp = -pgamma * rdamp * rdamp * rdamp;
                if (damp > -50.0) {
                  double expdamp = exp(damp);
                  scale3 = 1.0 - expdamp;
                  scale5 = 1.0 - expdamp * (1.0 - damp);
                  scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
                }

                final double dsc3 = scale3;
                final double dsc5 = scale5;
                final double dsc7 = scale7;
                final double rr3 = rr1 * rr2;
                final double rr5 = 3.0 * rr3 * rr2;
                final double rr7 = 5.0 * rr5 * rr2;
                final double drr3 = (1.0 - dsc3) * rr3;
                final double drr5 = (1.0 - dsc5) * rr5;
                final double drr7 = (1.0 - dsc7) * rr7;

                final double dkr = dkx * xr + dky * yr + dkz * zr;
                final double qkx = 2.0 * (qkxx * xr + qkxy * yr + qkxz * zr);
                final double qky = 2.0 * (qkxy * xr + qkyy * yr + qkyz * zr);
                final double qkz = 2.0 * (qkxz * xr + qkyz * yr + qkzz * zr);
                final double qkr = (qkx * xr + qky * yr + qkz * zr) * 0.5;
                final double bn123k = bn1 * ck - bn2 * dkr + bn3 * qkr;
                final double drr357k = drr3 * ck - drr5 * dkr + drr7 * qkr;
                final double fimx = -xr * bn123k - bn1 * dkx + bn2 * qkx;
                final double fimy = -yr * bn123k - bn1 * dky + bn2 * qky;
                final double fimz = -zr * bn123k - bn1 * dkz + bn2 * qkz;
                final double fidx = -xr * drr357k - drr3 * dkx + drr5 * qkx;
                final double fidy = -yr * drr357k - drr3 * dky + drr5 * qky;
                final double fidz = -zr * drr357k - drr3 * dkz + drr5 * qkz;

                final double dir = dix * xr + diy * yr + diz * zr;
                final double qix = 2.0 * (qixx * xr + qixy * yr + qixz * zr);
                final double qiy = 2.0 * (qixy * xr + qiyy * yr + qiyz * zr);
                final double qiz = 2.0 * (qixz * xr + qiyz * yr + qizz * zr);
                final double qir = (qix * xr + qiy * yr + qiz * zr) * 0.5;
                final double bn123i = bn1 * ci + bn2 * dir + bn3 * qir;
                final double ddr357i = drr3 * ci + drr5 * dir + drr7 * qir;
                final double fkmx = xr * bn123i - bn1 * dix - bn2 * qix;
                final double fkmy = yr * bn123i - bn1 * diy - bn2 * qiy;
                final double fkmz = zr * bn123i - bn1 * diz - bn2 * qiz;
                final double fkdx = xr * ddr357i - drr3 * dix - drr5 * qix;
                final double fkdy = yr * ddr357i - drr3 * diy - drr5 * qiy;
                final double fkdz = zr * ddr357i - drr3 * diz - drr5 * qiz;
                fix += selfScale * (fimx - fidx);
                fiy += selfScale * (fimy - fidy);
                fiz += selfScale * (fimz - fidz);
                final double xc = selfScale * (fkmx - fkdx);
                final double yc = selfScale * (fkmy - fkdy);
                final double zc = selfScale * (fkmz - fkdz);
                final double fkx = xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0];
                final double fky = xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1];
                final double fkz = xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2];
                field.add(threadID, k, fkx, fky, fkz);
                fieldCR.add(threadID, k, fkx, fky, fkz);
              }
            }
            field.add(threadID, i, fix, fiy, fiz);
            fieldCR.add(threadID, i, fix, fiy, fiz);
          }
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return permanentSchedule;
      }

      @Override
      public void start() {
        time = -System.nanoTime();
        threadID = getThreadIndex();
        count = 0;
        int nAtoms = atoms.length;
        if (inductionMaskLocal == null || inductionMaskLocal.length < nAtoms) {
          inductionMaskLocal = new double[nAtoms];
          energyMaskLocal = new double[nAtoms];
          fill(inductionMaskLocal, 1.0);
          fill(energyMaskLocal, 1.0);
        }
      }
    }

    private class PermanentRealSpaceFieldOctreeLoop extends IntegerForLoop {

      private final double[] dx_local;
      private final double[][] transOp;
      private int threadID;
      private double[] inductionMaskLocal;
      private double[] energyMaskLocal;
      private int count;
      private static final double theta = 0.5;

      PermanentRealSpaceFieldOctreeLoop() {
        super();
        dx_local = new double[3];
        transOp = new double[3][3];
      }

      public double[] getMomentsGeometric(Atom[] activeAtoms, NeighborList.Cell cell) {
        logger.info("** Enters getMomentsGeometric in PermFieldRegion");
        // Zero out total charge, dipole and quadrupole components.
        var netchg = 0.0;
        var netdpl = 0.0;
        var xdpl = 0.0;
        var ydpl = 0.0;
        var zdpl = 0.0;
        var xxqdp = 0.0;
        var xyqdp = 0.0;
        var xzqdp = 0.0;
        var yxqdp = 0.0;
        var yyqdp = 0.0;
        var yzqdp = 0.0;
        var zxqdp = 0.0;
        var zyqdp = 0.0;
        var zzqdp = 0.0;

        // Find the geometric center of the cell.
        int aind = cell.getIndices()[0];
        int bind = cell.getIndices()[1];
        int cind = cell.getIndices()[2];
        double xmid = cell.getSideLength() * (aind) + cell.getSideLength() * 0.5;
        double ymid = cell.getSideLength() * (bind) + cell.getSideLength() * 0.5;
        double zmid = cell.getSideLength() * (cind) + cell.getSideLength() * 0.5;

//    logger.info(format("Geometric Center of Cell %d %d %d : %4.4f , %4.4f , %4.4f",aind,bind,cind,xmid,ymid,zmid));

        // Get atom indices within Cell
        List<NeighborList.AtomIndex> atomIndexList = cell.getAtomIndexList();
        int n = atomIndexList.size();
        Atom[] cellAtoms = new Atom[n];
//    logger.info(format("cMG Number of atoms in cell %d %d %d : %d",aind,bind,cind,cellAtoms.length));

        // Get array of atoms contained in current cell
        int j = 0;
        for (NeighborList.AtomIndex index : atomIndexList){
          cellAtoms[j] = activeAtoms[index.i];
          j++;
        }

        double[] xcm = new double[n];
        double[] ycm = new double[n];
        double[] zcm = new double[n];
        int k = 0;
        for (Atom atom : cellAtoms) {
          xcm[k] = atom.getX() - xmid;
          ycm[k] = atom.getY() - ymid;
          zcm[k] = atom.getZ() - zmid;
          k++;
        }

        // Account for charge, dipoles and induced dipoles.
        k = 0;
        for (Atom atom : cellAtoms) {
          int i = atom.getIndex() - 1;
          double[] globalMultipolei = globalMultipole[0][i];
//      double[] inducedDipolei = inducedDipole[0][i];

          var ci = globalMultipolei[t000];
          var dix = globalMultipolei[t100];
          var diy = globalMultipolei[t010];
          var diz = globalMultipolei[t001];
//      var uix = inducedDipolei[0];
//      var uiy = inducedDipolei[1];
//      var uiz = inducedDipolei[2];

          netchg += ci;
          xdpl += xcm[k] * ci + dix;
          ydpl += ycm[k] * ci + diy;
          zdpl += zcm[k] * ci + diz;
          xxqdp += xcm[k] * xcm[k] * ci + 2.0 * xcm[k] * (dix);
          xyqdp += xcm[k] * ycm[k] * ci + xcm[k] * (diy) + ycm[k] * (dix);
          xzqdp += xcm[k] * zcm[k] * ci + xcm[k] * (diz) + zcm[k] * (dix);
          yxqdp += ycm[k] * xcm[k] * ci + ycm[k] * (dix) + xcm[k] * (diy);
          yyqdp += ycm[k] * ycm[k] * ci + 2.0 * ycm[k] * (diy);
          yzqdp += ycm[k] * zcm[k] * ci + ycm[k] * (diz) + zcm[k] * (diy);
          zxqdp += zcm[k] * xcm[k] * ci + zcm[k] * (dix) + xcm[k] * (diz);
          zyqdp += zcm[k] * ycm[k] * ci + zcm[k] * (diy) + ycm[k] * (diz);
          zzqdp += zcm[k] * zcm[k] * ci + 2.0 * zcm[k] * (diz);

//      logger.info(format("Cell %d %d %d    Charge %4.4f    Dipole %4.4f %4.4f %4.4f",aind,bind,cind,netchg,xdpl,ydpl,zdpl));

//      xdpl += xcm[k] * ci + dix + uix;
//      ydpl += ycm[k] * ci + diy + uiy;
//      zdpl += zcm[k] * ci + diz + uiz;
//      xxqdp += xcm[k] * xcm[k] * ci + 2.0 * xcm[k] * (dix + uix);
//      xyqdp += xcm[k] * ycm[k] * ci + xcm[k] * (diy + uiy) + ycm[k] * (dix + uix);
//      xzqdp += xcm[k] * zcm[k] * ci + xcm[k] * (diz + uiz) + zcm[k] * (dix + uix);
//      yxqdp += ycm[k] * xcm[k] * ci + ycm[k] * (dix + uix) + xcm[k] * (diy + uiy);
//      yyqdp += ycm[k] * ycm[k] * ci + 2.0 * ycm[k] * (diy + uiy);
//      yzqdp += ycm[k] * zcm[k] * ci + ycm[k] * (diz + uiz) + zcm[k] * (diy + uiy);
//      zxqdp += zcm[k] * xcm[k] * ci + zcm[k] * (dix + uix) + xcm[k] * (diz + uiz);
//      zyqdp += zcm[k] * ycm[k] * ci + zcm[k] * (diy + uiy) + ycm[k] * (diz + uiz);
//      zzqdp += zcm[k] * zcm[k] * ci + 2.0 * zcm[k] * (diz + uiz);

          k++;
        }

        // Convert the quadrupole from traced to traceless form.
        var qave = (xxqdp + yyqdp + zzqdp) / 3.0;
        xxqdp = 1.5 * (xxqdp - qave);
        xyqdp = 1.5 * xyqdp;
        xzqdp = 1.5 * xzqdp;
        yxqdp = 1.5 * yxqdp;
        yyqdp = 1.5 * (yyqdp - qave);
        yzqdp = 1.5 * yzqdp;
        zxqdp = 1.5 * zxqdp;
        zyqdp = 1.5 * zyqdp;
        zzqdp = 1.5 * (zzqdp - qave);

        // Add the traceless atomic quadrupoles to total quadrupole.
        for (Atom atom : cellAtoms) {
          int i = atom.getIndex() - 1;
          double[] globalMultipolei = globalMultipole[0][i];
          var qixx = globalMultipolei[t200];
          var qiyy = globalMultipolei[t020];
          var qizz = globalMultipolei[t002];
          var qixy = globalMultipolei[t110];
          var qixz = globalMultipolei[t101];
          var qiyz = globalMultipolei[t011];
          xxqdp += qixx;
          xyqdp += qixy;
          xzqdp += qixz;
          yxqdp += qixy;
          yyqdp += qiyy;
          yzqdp += qiyz;
          zxqdp += qixz;
          zyqdp += qiyz;
          zzqdp += qizz;
        }

        // Convert dipole to Debye and quadrupole to Buckingham.
        xdpl = xdpl * ELEC_ANG_TO_DEBYE;
        ydpl = ydpl * ELEC_ANG_TO_DEBYE;
        zdpl = zdpl * ELEC_ANG_TO_DEBYE;
        xxqdp = xxqdp * ELEC_ANG_TO_DEBYE;
        xyqdp = xyqdp * ELEC_ANG_TO_DEBYE;
        xzqdp = xzqdp * ELEC_ANG_TO_DEBYE;
        yxqdp = yxqdp * ELEC_ANG_TO_DEBYE;
        yyqdp = yyqdp * ELEC_ANG_TO_DEBYE;
        yzqdp = yzqdp * ELEC_ANG_TO_DEBYE;
        zxqdp = zxqdp * ELEC_ANG_TO_DEBYE;
        zyqdp = zyqdp * ELEC_ANG_TO_DEBYE;
        zzqdp = zzqdp * ELEC_ANG_TO_DEBYE;

        // Get dipole magnitude and diagonalize quadrupole tensor.
        netdpl = sqrt(xdpl * xdpl + ydpl * ydpl + zdpl * zdpl);
        double[][] a = new double[3][3];
        a[0][0] = xxqdp;
        a[0][1] = xyqdp;
        a[0][2] = xzqdp;
        a[1][0] = yxqdp;
        a[1][1] = yyqdp;
        a[1][2] = yzqdp;
        a[2][0] = zxqdp;
        a[2][1] = zyqdp;
        a[2][2] = zzqdp;
        EigenDecomposition e = new EigenDecomposition(new Array2DRowRealMatrix(a));
        // Eigenvalues are returned in descending order, but logged below in ascending order.
        var netqdp = e.getRealEigenvalues();

//    logger.info(format("\n Electric Moments for Cell A%d B%d C%d\n",cell.getIndices()[0], cell.getIndices()[1], cell.getIndices()[2]));
//    logger.info(format("  Total Electric Charge:    %13.5f Electrons\n", netchg));
//    logger.info(format("  Dipole Moment Magnitude:  %13.5f Debye\n", netdpl));
//    logger.info(format("  Dipole X,Y,Z-Components:  %13.5f %13.5f %13.5f\n", xdpl, ydpl, zdpl));
//    logger.info(format("  Quadrupole Moment Tensor: %13.5f %13.5f %13.5f", xxqdp, xyqdp, xzqdp));
//    logger.info(format("       (Buckinghams)        %13.5f %13.5f %13.5f", yxqdp, yyqdp, yzqdp));
//    logger.info(format("                            %13.5f %13.5f %13.5f\n", zxqdp, zyqdp, zzqdp));
//    logger.info(
//            format(
//                    "  Principal Axes Quadrupole %13.5f %13.5f %13.5f\n", netqdp[2], netqdp[1], netqdp[0]));

        return new double[]{xmid,ymid,zmid,netchg,xdpl,ydpl,zdpl,xxqdp,yyqdp,zzqdp,xyqdp,xzqdp,yzqdp};
      }

      @Override
      public void finish() {
        sharedCount.addAndGet(count);
      }

      @Override
      public void run(int lb, int ub) {
        int[][] lists = neighborLists[0];
        int[][] ewalds = realSpaceLists[0];
        int[] counts = realSpaceCounts[0];
        int[][] preLists = preconditionerLists[0];
        int[] preCounts = preconditionerCounts[0];
        final double[] x = coordinates[0][0];
        final double[] y = coordinates[0][1];
        final double[] z = coordinates[0][2];
        final double[][] mpole = globalMultipole[0];
        // Loop over atom chunk.
        for (int i = lb; i <= ub; i++) {
          if (!use[i]) {
            continue;
          }

          // Zero out accumulation variables for the field at atom i.
          double fix = 0.0;
          double fiy = 0.0;
          double fiz = 0.0;
          double fixCR = 0.0;
          double fiyCR = 0.0;
          double fizCR = 0.0;

          final int moleculei = molecule[i];
          final double pdi = ipdamp[i];
          final double pti = thole[i];
          final double xi = x[i];
          final double yi = y[i];
          final double zi = z[i];
          final double[] globalMultipolei = mpole[i];
          final double ci = globalMultipolei[0];
          final double dix = globalMultipolei[t100];
          final double diy = globalMultipolei[t010];
          final double diz = globalMultipolei[t001];
          final double qixx = globalMultipolei[t200] * oneThird;
          final double qiyy = globalMultipolei[t020] * oneThird;
          final double qizz = globalMultipolei[t002] * oneThird;
          final double qixy = globalMultipolei[t110] * oneThird;
          final double qixz = globalMultipolei[t101] * oneThird;
          final double qiyz = globalMultipolei[t011] * oneThird;

          // Apply field masking rules.
          applyMask(i, null, energyMaskLocal, inductionMaskLocal);

          // Loop over the neighbor list.
          int[] list = lists[i];
          counts[i] = 0;
          preCounts[i] = 0;
          int[] ewald = ewalds[i];
          int[] preList = preLists[i];

          boolean cellMultipole = false;
          // Loop over cells
          for (int l = 0; l < cells.length; l++){
            for (int m = 0; m < cells.length; m++){
              for (int n = 0; n < cells.length; n++){
                NeighborList.Cell currentCell = cells[l][m][n];
                double[] cellComponents = getMomentsGeometric(atoms,currentCell);
                double xc = cellComponents[0];
                double yc = cellComponents[1];
                double zc = cellComponents[2];

                double distance = Math.sqrt(Math.pow(xc - xi,2) + Math.pow(yc - yi,2) + Math.pow(zc - zi,2));
                if (currentCell.getSideLength()/distance < theta){
                  // Cell is far enough away to consider a single cell-based multipole interaction
                  list = new int[]{0};
                  cellMultipole = true;
                } else {
                  // Cell is too close for cell-based multipole interactions - compute individual atomic multipole interactions
                  list = currentCell.getAtomIndices();
                  cellMultipole = false;
                }

                for (int k : list) {
                  if (!use[k] && !cellMultipole) {
                    continue;
                  }
                  if (lambdaMode == LambdaMode.VAPOR) {
                    boolean sameMolecule = (moleculei == molecule[k]);
                    if ((intermolecularSoftcore && !sameMolecule)
                            || (intramolecularSoftcore && sameMolecule)) {
                      continue;
                    }
                  }

                  double tmpx;
                  double tmpy;
                  double tmpz;
                  if(cellMultipole){
                    tmpx = xc;
                    tmpy = yc;
                    tmpz = zc;
                  } else {
                    tmpx = x[k];
                    tmpy = y[k];
                    tmpz = z[k];
                  }

                  final double xk = tmpx;
                  final double yk = tmpy;
                  final double zk = tmpz;
                  dx_local[0] = xk - xi;
                  dx_local[1] = yk - yi;
                  dx_local[2] = zk - zi;
                  final double r2 = crystal.image(dx_local);
                  if (r2 <= off2) {
                    count++;
                    // Store a short neighbor list for the SCF.
                    if (ewald.length <= counts[i]) {
                      int len = ewald.length;
                      ewalds[i] = copyOf(ewald, len + 10);
                      ewald = ewalds[i];
                    }
                    // TODO: check what ewald array is used for and how k is referenced
                    ewald[counts[i]++] = k;
                    final double xr = dx_local[0];
                    final double yr = dx_local[1];
                    final double zr = dx_local[2];
                    final double pdk = ipdamp[k];
                    final double ptk = thole[k];
                    final double[] globalMultipolek = mpole[k];
                    final double ck = globalMultipolek[t000];
                    final double dkx = globalMultipolek[t100];
                    final double dky = globalMultipolek[t010];
                    final double dkz = globalMultipolek[t001];
                    final double qkxx = globalMultipolek[t200] * oneThird;
                    final double qkyy = globalMultipolek[t020] * oneThird;
                    final double qkzz = globalMultipolek[t002] * oneThird;
                    final double qkxy = globalMultipolek[t110] * oneThird;
                    final double qkxz = globalMultipolek[t101] * oneThird;
                    final double qkyz = globalMultipolek[t011] * oneThird;
                    double r = sqrt(r2);

                    // Store a short neighbor list for the SCF pre-conditioner.
                    if (r < preconditionerCutoff) {
                      if (preList.length <= preCounts[i]) {
                        int len = preList.length;
                        preLists[i] = copyOf(preList, len + 10);
                        preList = preLists[i];
                      }
                      // TODO: check what preList array is used for and how k is referenced
                      preList[preCounts[i]++] = k;
                    }

                    // Calculate the error function damping terms.
                    final double ralpha = aewald * r;
                    final double exp2a = exp(-ralpha * ralpha);
                    final double rr1 = 1.0 / r;
                    final double rr2 = rr1 * rr1;
                    final double bn0 = erfc(ralpha) * rr1;
                    final double bn1 = (bn0 + an0 * exp2a) * rr2;
                    final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                    final double bn3 = (5.0 * bn2 + an2 * exp2a) * rr2;

                    // Compute the error function scaled and unscaled terms.
                    double scale3 = 1.0;
                    double scale5 = 1.0;
                    double scale7 = 1.0;
                    double damp = pdi * pdk;
                    final double pgamma = min(pti, ptk);
                    final double rdamp = r * damp;
                    damp = -pgamma * rdamp * rdamp * rdamp;
                    if (damp > -50.0) {
                      double expdamp = exp(damp);
                      scale3 = 1.0 - expdamp;
                      scale5 = 1.0 - expdamp * (1.0 - damp);
                      scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
                    }
                    final double scale = inductionMaskLocal[k];
                    final double scalep = energyMaskLocal[k];
                    // Thole damping multiplied by the group-based mask.
                    final double dsc3 = scale3 * scale;
                    final double dsc5 = scale5 * scale;
                    final double dsc7 = scale7 * scale;
                    // Thole damping multiplied by the energy mask.
                    final double psc3 = scale3 * scalep;
                    final double psc5 = scale5 * scalep;
                    final double psc7 = scale7 * scalep;
                    final double rr3 = rr1 * rr2;
                    final double rr5 = 3.0 * rr3 * rr2;
                    final double rr7 = 5.0 * rr5 * rr2;
                    // 1.0 minus induction masks and Thole damping.
                    final double drr3 = (1.0 - dsc3) * rr3;
                    final double drr5 = (1.0 - dsc5) * rr5;
                    final double drr7 = (1.0 - dsc7) * rr7;
                    // 1.0 minus energy masks and Thole damping.
                    final double prr3 = (1.0 - psc3) * rr3;
                    final double prr5 = (1.0 - psc5) * rr5;
                    final double prr7 = (1.0 - psc7) * rr7;
                    final double dir = dix * xr + diy * yr + diz * zr;
                    final double qix = 2.0 * (qixx * xr + qixy * yr + qixz * zr);
                    final double qiy = 2.0 * (qixy * xr + qiyy * yr + qiyz * zr);
                    final double qiz = 2.0 * (qixz * xr + qiyz * yr + qizz * zr);
                    final double qir = (qix * xr + qiy * yr + qiz * zr) * 0.5;
                    // Ewald field for atom k (no masking).
                    final double bn123i = bn1 * ci + bn2 * dir + bn3 * qir;
                    final double fkmx = xr * bn123i - bn1 * dix - bn2 * qix;
                    final double fkmy = yr * bn123i - bn1 * diy - bn2 * qiy;
                    final double fkmz = zr * bn123i - bn1 * diz - bn2 * qiz;
                    // Correct Ewald field for over-counted induction interactions.
                    final double ddr357i = drr3 * ci + drr5 * dir + drr7 * qir;
                    final double fkdx = xr * ddr357i - drr3 * dix - drr5 * qix;
                    final double fkdy = yr * ddr357i - drr3 * diy - drr5 * qiy;
                    final double fkdz = zr * ddr357i - drr3 * diz - drr5 * qiz;
                    field.add(threadID, k, fkmx - fkdx, fkmy - fkdy, fkmz - fkdz);
                    // Correct Ewald field for over-counted energy interactions.
                    final double prr357i = prr3 * ci + prr5 * dir + prr7 * qir;
                    final double fkpx = xr * prr357i - prr3 * dix - prr5 * qix;
                    final double fkpy = yr * prr357i - prr3 * diy - prr5 * qiy;
                    final double fkpz = zr * prr357i - prr3 * diz - prr5 * qiz;
                    fieldCR.add(threadID, k, fkmx - fkpx, fkmy - fkpy, fkmz - fkpz);
                    final double dkr = dkx * xr + dky * yr + dkz * zr;
                    final double qkx = 2.0 * (qkxx * xr + qkxy * yr + qkxz * zr);
                    final double qky = 2.0 * (qkxy * xr + qkyy * yr + qkyz * zr);
                    final double qkz = 2.0 * (qkxz * xr + qkyz * yr + qkzz * zr);
                    final double qkr = (qkx * xr + qky * yr + qkz * zr) * 0.5;
                    final double bn123k = bn1 * ck - bn2 * dkr + bn3 * qkr;
                    // Ewald field for atom i (no masking).
                    final double fimx = -xr * bn123k - bn1 * dkx + bn2 * qkx;
                    final double fimy = -yr * bn123k - bn1 * dky + bn2 * qky;
                    final double fimz = -zr * bn123k - bn1 * dkz + bn2 * qkz;
                    // Correct Ewald field for over-counted induction interactions.
                    final double drr357k = drr3 * ck - drr5 * dkr + drr7 * qkr;
                    final double fidx = -xr * drr357k - drr3 * dkx + drr5 * qkx;
                    final double fidy = -yr * drr357k - drr3 * dky + drr5 * qky;
                    final double fidz = -zr * drr357k - drr3 * dkz + drr5 * qkz;
                    fix += fimx - fidx;
                    fiy += fimy - fidy;
                    fiz += fimz - fidz;
                    // Correct Ewald field for over-counted energy interactions.
                    final double prr357k = prr3 * ck - prr5 * dkr + prr7 * qkr;
                    final double fipx = -xr * prr357k - prr3 * dkx + prr5 * qkx;
                    final double fipy = -yr * prr357k - prr3 * dky + prr5 * qky;
                    final double fipz = -zr * prr357k - prr3 * dkz + prr5 * qkz;
                    fixCR += fimx - fipx;
                    fiyCR += fimy - fipy;
                    fizCR += fimz - fipz;
                  }
                }

              }
            }
          }


          // Add in field contributions at Atom i.
          field.add(threadID, i, fix, fiy, fiz);
          fieldCR.add(threadID, i, fixCR, fiyCR, fizCR);
          // Remove field masking rules.
          removeMask(i, null, energyMaskLocal, inductionMaskLocal);
        }

        // Loop over symmetry mates.
        List<SymOp> symOps = crystal.spaceGroup.symOps;
        int nSymm = symOps.size();
        for (int iSymm = 1; iSymm < nSymm; iSymm++) {
          SymOp symOp = crystal.spaceGroup.getSymOp(iSymm);
          crystal.getTransformationOperator(symOp, transOp);
          lists = neighborLists[iSymm];
          ewalds = realSpaceLists[iSymm];
          counts = realSpaceCounts[iSymm];
          preLists = preconditionerLists[iSymm];
          preCounts = preconditionerCounts[iSymm];
          double[] xs = coordinates[iSymm][0];
          double[] ys = coordinates[iSymm][1];
          double[] zs = coordinates[iSymm][2];
          double[][] mpoles = globalMultipole[iSymm];

          // Loop over atoms in a chunk of the asymmetric unit.
          for (int i = lb; i <= ub; i++) {
            if (!use[i]) {
              continue;
            }
            // Zero out accumulation variables for the field at atom i.
            double fix = 0.0;
            double fiy = 0.0;
            double fiz = 0.0;
            final double pdi = ipdamp[i];
            final double pti = thole[i];
            final double[] multipolei = mpole[i];
            final double ci = multipolei[t000];
            final double dix = multipolei[t100];
            final double diy = multipolei[t010];
            final double diz = multipolei[t001];
            final double qixx = multipolei[t200] * oneThird;
            final double qiyy = multipolei[t020] * oneThird;
            final double qizz = multipolei[t002] * oneThird;
            final double qixy = multipolei[t110] * oneThird;
            final double qixz = multipolei[t101] * oneThird;
            final double qiyz = multipolei[t011] * oneThird;
            final double xi = x[i];
            final double yi = y[i];
            final double zi = z[i];

            // Loop over the neighbor list.
            final int[] list = lists[i];
            counts[i] = 0;
            preCounts[i] = 0;
            int[] ewald = ewalds[i];
            int[] preList = preLists[i];
            for (int k : list) {
              if (!use[k]) {
                continue;
              }
              final double xk = xs[k];
              final double yk = ys[k];
              final double zk = zs[k];
              dx_local[0] = xk - xi;
              dx_local[1] = yk - yi;
              dx_local[2] = zk - zi;
              final double r2 = crystal.image(dx_local);
              if (r2 <= off2) {
                count++;
                // Store a short neighbor list for the SCF.
                if (ewald.length <= counts[i]) {
                  int len = ewald.length;
                  ewalds[i] = copyOf(ewald, len + 10);
                  ewald = ewalds[i];
                }
                ewald[counts[i]++] = k;
                double selfScale = 1.0;
                if (i == k) {
                  selfScale = 0.5;
                }
                final double xr = dx_local[0];
                final double yr = dx_local[1];
                final double zr = dx_local[2];
                final double pdk = ipdamp[k];
                final double ptk = thole[k];
                final double[] multipolek = mpoles[k];
                final double ck = multipolek[t000];
                final double dkx = multipolek[t100];
                final double dky = multipolek[t010];
                final double dkz = multipolek[t001];
                final double qkxx = multipolek[t200] * oneThird;
                final double qkyy = multipolek[t020] * oneThird;
                final double qkzz = multipolek[t002] * oneThird;
                final double qkxy = multipolek[t110] * oneThird;
                final double qkxz = multipolek[t101] * oneThird;
                final double qkyz = multipolek[t011] * oneThird;
                final double r = sqrt(r2);
                if (r < preconditionerCutoff) {
                  if (preList.length <= preCounts[i]) {
                    int len = preList.length;
                    preLists[i] = copyOf(preList, len + 10);
                    preList = preLists[i];
                  }
                  preList[preCounts[i]++] = k;
                }

                // Calculate the error function damping terms.
                final double ralpha = aewald * r;
                final double exp2a = exp(-ralpha * ralpha);
                final double rr1 = 1.0 / r;
                final double rr2 = rr1 * rr1;
                final double bn0 = erfc(ralpha) * rr1;
                final double bn1 = (bn0 + an0 * exp2a) * rr2;
                final double bn2 = (3.0 * bn1 + an1 * exp2a) * rr2;
                final double bn3 = (5.0 * bn2 + an2 * exp2a) * rr2;

                // Compute the error function scaled and unscaled terms.
                double scale3 = 1.0;
                double scale5 = 1.0;
                double scale7 = 1.0;
                double damp = pdi * pdk;
                final double pgamma = min(pti, ptk);
                final double rdamp = r * damp;
                damp = -pgamma * rdamp * rdamp * rdamp;
                if (damp > -50.0) {
                  double expdamp = exp(damp);
                  scale3 = 1.0 - expdamp;
                  scale5 = 1.0 - expdamp * (1.0 - damp);
                  scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
                }

                final double dsc3 = scale3;
                final double dsc5 = scale5;
                final double dsc7 = scale7;
                final double rr3 = rr1 * rr2;
                final double rr5 = 3.0 * rr3 * rr2;
                final double rr7 = 5.0 * rr5 * rr2;
                final double drr3 = (1.0 - dsc3) * rr3;
                final double drr5 = (1.0 - dsc5) * rr5;
                final double drr7 = (1.0 - dsc7) * rr7;

                final double dkr = dkx * xr + dky * yr + dkz * zr;
                final double qkx = 2.0 * (qkxx * xr + qkxy * yr + qkxz * zr);
                final double qky = 2.0 * (qkxy * xr + qkyy * yr + qkyz * zr);
                final double qkz = 2.0 * (qkxz * xr + qkyz * yr + qkzz * zr);
                final double qkr = (qkx * xr + qky * yr + qkz * zr) * 0.5;
                final double bn123k = bn1 * ck - bn2 * dkr + bn3 * qkr;
                final double drr357k = drr3 * ck - drr5 * dkr + drr7 * qkr;
                final double fimx = -xr * bn123k - bn1 * dkx + bn2 * qkx;
                final double fimy = -yr * bn123k - bn1 * dky + bn2 * qky;
                final double fimz = -zr * bn123k - bn1 * dkz + bn2 * qkz;
                final double fidx = -xr * drr357k - drr3 * dkx + drr5 * qkx;
                final double fidy = -yr * drr357k - drr3 * dky + drr5 * qky;
                final double fidz = -zr * drr357k - drr3 * dkz + drr5 * qkz;

                final double dir = dix * xr + diy * yr + diz * zr;
                final double qix = 2.0 * (qixx * xr + qixy * yr + qixz * zr);
                final double qiy = 2.0 * (qixy * xr + qiyy * yr + qiyz * zr);
                final double qiz = 2.0 * (qixz * xr + qiyz * yr + qizz * zr);
                final double qir = (qix * xr + qiy * yr + qiz * zr) * 0.5;
                final double bn123i = bn1 * ci + bn2 * dir + bn3 * qir;
                final double ddr357i = drr3 * ci + drr5 * dir + drr7 * qir;
                final double fkmx = xr * bn123i - bn1 * dix - bn2 * qix;
                final double fkmy = yr * bn123i - bn1 * diy - bn2 * qiy;
                final double fkmz = zr * bn123i - bn1 * diz - bn2 * qiz;
                final double fkdx = xr * ddr357i - drr3 * dix - drr5 * qix;
                final double fkdy = yr * ddr357i - drr3 * diy - drr5 * qiy;
                final double fkdz = zr * ddr357i - drr3 * diz - drr5 * qiz;
                fix += selfScale * (fimx - fidx);
                fiy += selfScale * (fimy - fidy);
                fiz += selfScale * (fimz - fidz);
                final double xc = selfScale * (fkmx - fkdx);
                final double yc = selfScale * (fkmy - fkdy);
                final double zc = selfScale * (fkmz - fkdz);
                final double fkx = xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0];
                final double fky = xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1];
                final double fkz = xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2];
                field.add(threadID, k, fkx, fky, fkz);
                fieldCR.add(threadID, k, fkx, fky, fkz);
              }
            }
            field.add(threadID, i, fix, fiy, fiz);
            fieldCR.add(threadID, i, fix, fiy, fiz);
          }
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return permanentSchedule;
      }

      @Override
      public void start() {
        threadID = getThreadIndex();
        count = 0;
        int nAtoms = atoms.length;
        if (inductionMaskLocal == null || inductionMaskLocal.length < nAtoms) {
          inductionMaskLocal = new double[nAtoms];
          energyMaskLocal = new double[nAtoms];
          fill(inductionMaskLocal, 1.0);
          fill(energyMaskLocal, 1.0);
        }
      }
    }

  }
}
