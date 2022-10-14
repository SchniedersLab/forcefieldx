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
package ffx.potential.nonbonded.pme;

import static ffx.numerics.special.Erf.erfc;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwald.EwaldParameters;
import ffx.potential.parameters.ForceField;
import ffx.potential.utils.EnergyException;
import ffx.utilities.Constants;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Parallel pre-conditioned conjugate gradient solver for the self-consistent field.
 * <p>
 * This solves the linear system A x = b using an iterative approach where A: is an N x N symmetric,
 * positive definite matrix B: is a known vector of length N. x: is the unknown vector of length N.
 * <p>
 * For the AMOEBA SCF, the linear system Ax = b is usually denoted: C u = E_multipoles where C = [
 * alpha^-1 - T ] u are the induced dipoles. E_multipoles is the direct field from permanent
 * multipoles.
 * <p>
 * The matrix alpha^-1 is the inverse of the N x N diagonal polarizability matrix. The matrix T is
 * the N x N matrix that produces the field due to induced dipoles.
 * <p>
 * Initialization: 1) Compute the residual:        r_0 = E_mutipoles - C u_direct = E_direct_induced
 * 2) Compute the preconditioner:  z_0 = T r_0 3) Compute the conjugate:       p_0 = z_0 4) Initial
 * loop index:          k = 0
 * <p>
 * Then loop over: 1) Update the step size:        alpha = r_k dot z_k / (p_k C p_k) 2) Update the
 * induced dipoles:  u_k+1 = u_k + alpha p_k 3) Update the residual:         r_k+1 = r_k - alpha C
 * p_k 4) Check for convergence of r_k+1. 5) Update the preconditioner:   z_k+1 = T r_k+1 6) Update
 * the step size:        beta = r_k+1 dot z_k+1 / (r_k dot z_k) 7) Update the conjugate:        p_k+1
 * = z_k+1 + beta p_k 8) Update the loop index:       k = k + 1
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PCGSolver {

  private static final Logger logger = Logger.getLogger(PCGSolver.class.getName());

  private final InitResidualRegion initResidualRegion;
  private final InitConjugateRegion initConjugateRegion;
  private final UpdateResidualRegion updateResidualRegion;
  private final UpdateConjugateRegion updateConjugateRegion;
  private final PreconditionerRegion preconditionerRegion;
  private final double poleps;
  public final double preconditionerCutoff;
  public final double preconditionerEwald;
  /**
   * Acceleration factor for induced dipole SCF iterations. This parameter weights the diagonol
   * elements of the preconditioning matrix.
   */
  public final double preconditionerScale;
  /**
   * Neighbor lists, without atoms beyond the preconditioner cutoff.
   * [nSymm][nAtoms][nIncludedNeighbors]
   */
  int[][][] preconditionerLists;
  /** Number of neighboring atoms within the preconditioner cutoff. [nSymm][nAtoms] */
  int[][] preconditionerCounts;
  /**
   * Residual vector (an electric field).
   */
  private double[][] r;
  /**
   * Residual vector for the chain-rule dipoles (an electric field).
   */
  private double[][] rCR;
  /**
   * Preconditioner dipoles (z = M^-1 r).
   */
  private double[][] z;
  /**
   * Preconditioner dipoles for the chain-rule dipoles (zCR = M^-1 rCR).
   */
  private double[][] zCR;
  /**
   * Conjugate search direction (induced dipoles).
   */
  private double[][] p;
  /**
   * Conjugate search direction for the chain-rule dipoles (induced dipoles).
   */
  private double[][] pCR;
  /** Work vector. */
  private double[][] vec;
  /** Work vector for the chain-rule dipoles. */
  private double[][] vecCR;
  /** An ordered array of atoms in the system. */
  private Atom[] atoms;
  /** Dimensions of [nsymm][xyz][nAtoms]. */
  private double[][][] coordinates;
  private double[] polarizability;
  private double[] ipdamp;
  private double[] thole;
  /**
   * When computing the polarization energy at Lambda there are 3 pieces.
   *
   * <p>1.) Upol(1) = The polarization energy computed normally (i.e. system with ligand).
   *
   * <p>2.) Uenv = The polarization energy of the system without the ligand.
   *
   * <p>3.) Uligand = The polarization energy of the ligand by itself.
   *
   * <p>Upol(L) = L*Upol(1) + (1-L)*(Uenv + Uligand)
   *
   * <p>Set the "use" array to true for all atoms for part 1.
   * <p>Set the "use" array to true for all atoms except the ligand for part 2.
   * <p>Set the "use" array to true only for the ligand atoms for part 3.
   *
   * <p>The "use" array can also be employed to turn off atoms for computing the electrostatic
   * energy of sub-structures.
   */
  private boolean[] use;
  /** Unit cell and spacegroup information. */
  private Crystal crystal;
  /** Dimensions of [nsymm][nAtoms][3] */
  private double[][][] inducedDipole;
  private double[][][] inducedDipoleCR;
  /** Direct induced dipoles. */
  private double[][] directDipole;
  private double[][] directDipoleCR;

  /** Field array. */
  private AtomicDoubleArray3D field;
  /** Chain rule field array. */
  private AtomicDoubleArray3D fieldCR;
  private EwaldParameters ewaldParameters;
  /**
   * The default ParallelTeam encapsulates the maximum number of threads used to parallelize the
   * electrostatics calculation.
   */
  private ParallelTeam parallelTeam;
  /** Pairwise schedule for load balancing. */
  private IntegerSchedule realSpaceSchedule;
  private long[] realSpaceSCFTime;

  /**
   * Constructor the PCG solver.
   *
   * @param maxThreads Number of threads.
   * @param poleps Convergence criteria (RMS Debye).
   * @param forceField Force field in use.
   * @param nAtoms Initial number of atoms.
   */
  public PCGSolver(int maxThreads, double poleps, ForceField forceField, int nAtoms) {
    this.poleps = poleps;
    preconditionerRegion = new PreconditionerRegion(maxThreads);
    initResidualRegion = new InitResidualRegion(maxThreads);
    initConjugateRegion = new InitConjugateRegion(maxThreads);
    updateResidualRegion = new UpdateResidualRegion(maxThreads);
    updateConjugateRegion = new UpdateConjugateRegion(maxThreads);

    // The size of the preconditioner neighbor list depends on the size of the preconditioner
    // cutoff.
    boolean preconditioner = forceField.getBoolean("USE_SCF_PRECONDITIONER", true);
    if (preconditioner) {
      preconditionerCutoff = forceField.getDouble("CG_PRECONDITIONER_CUTOFF", 4.5);
      preconditionerEwald = forceField.getDouble("CG_PRECONDITIONER_EWALD", 0.0);
      preconditionerScale = forceField.getDouble("CG_PRECONDITIONER_SCALE", 2.0);
    } else {
      preconditionerCutoff = 0.0;
      preconditionerEwald = 0.0;
      preconditionerScale = 0.0;
    }

    allocateVectors(nAtoms);
  }

  /**
   * Allocate storage for pre-conditioner neighbor list.
   *
   * @param nSymm Number of symmetry operators.
   * @param nAtoms Number of atoms.
   */
  public void allocateLists(int nSymm, int nAtoms) {
    int preconditionerListSize = 50;
    preconditionerLists = new int[nSymm][nAtoms][preconditionerListSize];
    preconditionerCounts = new int[nSymm][nAtoms];
  }

  /**
   * Allocate PCG vectors.
   *
   * @param nAtoms The number of atoms.
   */
  public void allocateVectors(int nAtoms) {
    if (r == null || r[0].length != nAtoms) {
      r = new double[3][nAtoms];
      rCR = new double[3][nAtoms];
      z = new double[3][nAtoms];
      zCR = new double[3][nAtoms];
      p = new double[3][nAtoms];
      pCR = new double[3][nAtoms];
      vec = new double[3][nAtoms];
      vecCR = new double[3][nAtoms];
    }
  }

  public void init(
      Atom[] atoms,
      double[][][] coordinates,
      double[] polarizability,
      double[] ipdamp,
      double[] thole,
      boolean[] use,
      Crystal crystal,
      double[][][] inducedDipole,
      double[][][] inducedDipoleCR,
      double[][] directDipole,
      double[][] directDipoleCR,
      AtomicDoubleArray3D field,
      AtomicDoubleArray3D fieldCR,
      EwaldParameters ewaldParameters,
      ParallelTeam parallelTeam,
      IntegerSchedule realSpaceSchedule,
      long[] realSpaceSCFTime) {
    this.atoms = atoms;
    this.coordinates = coordinates;
    this.polarizability = polarizability;
    this.ipdamp = ipdamp;
    this.thole = thole;
    this.use = use;
    this.crystal = crystal;
    this.inducedDipole = inducedDipole;
    this.inducedDipoleCR = inducedDipoleCR;
    this.directDipole = directDipole;
    this.directDipoleCR = directDipoleCR;
    this.field = field;
    this.fieldCR = fieldCR;
    this.ewaldParameters = ewaldParameters;
    this.parallelTeam = parallelTeam;
    this.realSpaceSchedule = realSpaceSchedule;
    this.realSpaceSCFTime = realSpaceSCFTime;
  }

  public int scfByPCG(boolean print, long startTime, ParticleMeshEwald pme) {
    long directTime = System.nanoTime() - startTime;
    // A request of 0 SCF cycles simplifies mutual polarization to direct polarization.
    StringBuilder sb = null;
    if (print) {
      sb = new StringBuilder("\n Self-Consistent Field\n Iter  RMS Change (Debye)  Time\n");
    }

    // Find the induced dipole field due to direct dipoles
    // (or predicted induced dipoles from previous steps).
    pme.computeInduceDipoleField();

    try {
      // Set the initial residual.
      parallelTeam.execute(initResidualRegion);

      // Compute the induced field due to the residual using short cut-offs.
      computePreconditioner(pme);

      // Set initial conjugate vector.
      parallelTeam.execute(initConjugateRegion);
    } catch (Exception e) {
      String message = "Exception initializing preconditioned CG.";
      logger.log(Level.SEVERE, message, e);
    }

    // Conjugate gradient iteration of the mutual induced dipoles.
    int completedSCFCycles = 0;
    int maxSCFCycles = 1000;
    double eps = 100.0;
    double previousEps;
    boolean done = false;
    while (!done) {
      long cycleTime = -System.nanoTime();

      // Store a copy of the current induced dipoles,
      // then set the induced dipoles to the conjugate vector.
      int nAtoms = atoms.length;
      for (int i = 0; i < nAtoms; i++) {
        if (use[i]) {
          vec[0][i] = inducedDipole[0][i][0];
          vec[1][i] = inducedDipole[0][i][1];
          vec[2][i] = inducedDipole[0][i][2];
          inducedDipole[0][i][0] = p[0][i];
          inducedDipole[0][i][1] = p[1][i];
          inducedDipole[0][i][2] = p[2][i];
          vecCR[0][i] = inducedDipoleCR[0][i][0];
          vecCR[1][i] = inducedDipoleCR[0][i][1];
          vecCR[2][i] = inducedDipoleCR[0][i][2];
          inducedDipoleCR[0][i][0] = pCR[0][i];
          inducedDipoleCR[0][i][1] = pCR[1][i];
          inducedDipoleCR[0][i][2] = pCR[2][i];
        } else {
          vec[0][i] = 0.0;
          vec[1][i] = 0.0;
          vec[2][i] = 0.0;
          inducedDipole[0][i][0] = 0.0;
          inducedDipole[0][i][1] = 0.0;
          inducedDipole[0][i][2] = 0.0;
          vecCR[0][i] = 0.0;
          vecCR[1][i] = 0.0;
          vecCR[2][i] = 0.0;
          inducedDipoleCR[0][i][0] = 0.0;
          inducedDipoleCR[0][i][1] = 0.0;
          inducedDipoleCR[0][i][2] = 0.0;
        }
      }

      // Find the induced dipole field due to the conjugate search direction p.
      pme.computeInduceDipoleField();

      try {
        /*
         * 1) Using the induced dipole field, compute Ap.
         * 2) Compute the step size alpha.
         * 3) Update the residual and induced dipole solution.
         */
        parallelTeam.execute(updateResidualRegion);

        // Compute the induced field due to the residual using short cut-offs.
        computePreconditioner(pme);

        /*
         * 1) Compute the dot product of the residual (r) and preconditioner field (z).
         * 2) Compute the step size beta.
         * 3) Update the conjugate search direction (p).
         */
        updateConjugateRegion.previousRDotZ = updateResidualRegion.getRDotZ();
        updateConjugateRegion.previousRDotZCR = updateResidualRegion.getRDotZCR();
        parallelTeam.execute(updateConjugateRegion);
      } catch (Exception e) {
        String message = "Exception during CG iteration.";
        logger.log(Level.SEVERE, message, e);
      }

      previousEps = eps;
      eps = max(updateConjugateRegion.getEps(), updateConjugateRegion.getEpsCR());
      completedSCFCycles++;
      eps = Constants.ELEC_ANG_TO_DEBYE * sqrt(eps / (double) nAtoms);
      cycleTime += System.nanoTime();
      if (print) {
        sb.append(
            format(
                " %4d     %15.10f %7.4f\n", completedSCFCycles, eps, cycleTime * Constants.NS2SEC));
      }

      // If the RMS Debye change increases, fail the SCF process.
      if (eps > previousEps) {
        if (sb != null) {
          logger.warning(sb.toString());
        }
        String message =
            format("Fatal SCF convergence failure: (%10.5f > %10.5f)\n", eps, previousEps);
        throw new EnergyException(message);
      }

      // The SCF should converge well before the max iteration check. Otherwise, fail the SCF
      // process.
      if (completedSCFCycles >= maxSCFCycles) {
        if (sb != null) {
          logger.warning(sb.toString());
        }
        String message = format("Maximum SCF iterations reached: (%d)\n", completedSCFCycles);
        throw new EnergyException(message);
      }

      // Check if the convergence criteria has been achieved.
      if (eps < poleps) {
        done = true;
      }
    }
    if (print) {
      sb.append(format(" Direct:                  %7.4f\n", Constants.NS2SEC * directTime));
      startTime = System.nanoTime() - startTime;
      sb.append(format(" Total:                   %7.4f", startTime * Constants.NS2SEC));
      logger.info(sb.toString());
    }

    // Find the final induced dipole field.
    pme.computeInduceDipoleField();

    return completedSCFCycles;
  }

  /**
   * Compute the induced field due to the Uind = polarizability * Eresidual.
   *
   * @param particleMeshEwald An instance of ParticleMeshEwald.
   */
  private void computePreconditioner(ParticleMeshEwald particleMeshEwald) {
    try {
      particleMeshEwald.expandInducedDipoles();
      // Use a special Ewald coefficient for the pre-conditioner.
      double aewaldTemp = ewaldParameters.aewald;
      ewaldParameters.setEwaldParameters(ewaldParameters.off, preconditionerEwald);
      int nAtoms = atoms.length;
      field.reset(parallelTeam, 0, nAtoms - 1);
      fieldCR.reset(parallelTeam, 0, nAtoms - 1);
      parallelTeam.execute(preconditionerRegion);
      field.reduce(parallelTeam, 0, nAtoms - 1);
      fieldCR.reduce(parallelTeam, 0, nAtoms - 1);
      ewaldParameters.setEwaldParameters(ewaldParameters.off, aewaldTemp);
    } catch (Exception e) {
      String message = "Exception computing the induced field for the preconditioner.";
      logger.log(Level.SEVERE, message, e);
    }
  }

  /**
   * Set the initial Residual r_0 = E_perm - u_0 / alpha + E_u_0.
   */
  private class InitResidualRegion extends ParallelRegion {

    private final InitResidualLoop[] initResidualLoops;

    public InitResidualRegion(int nt) {
      initResidualLoops = new InitResidualLoop[nt];
    }

    @Override
    public void run() throws Exception {
      try {
        int ti = getThreadIndex();
        if (initResidualLoops[ti] == null) {
          initResidualLoops[ti] = new InitResidualLoop();
        }
        int nAtoms = atoms.length;
        execute(0, nAtoms - 1, initResidualLoops[ti]);
      } catch (Exception e) {
        String message =
            "Fatal exception computing the mutual induced dipoles in thread "
                + getThreadIndex()
                + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }

    private class InitResidualLoop extends IntegerForLoop {

      @Override
      public void run(int lb, int ub) throws Exception {
        for (int i = lb; i <= ub; i++) {
          // Set initial residual to the direct field.
          if (use[i] && polarizability[i] > 0.0) {
            double ipolar = 1.0 / polarizability[i];
            r[0][i] = (directDipole[i][0] - inducedDipole[0][i][0]) * ipolar + field.getX(i);
            r[1][i] = (directDipole[i][1] - inducedDipole[0][i][1]) * ipolar + field.getY(i);
            r[2][i] = (directDipole[i][2] - inducedDipole[0][i][2]) * ipolar + field.getZ(i);
            rCR[0][i] = (directDipoleCR[i][0] - inducedDipoleCR[0][i][0]) * ipolar + fieldCR.getX(i);
            rCR[1][i] = (directDipoleCR[i][1] - inducedDipoleCR[0][i][1]) * ipolar + fieldCR.getY(i);
            rCR[2][i] = (directDipoleCR[i][2] - inducedDipoleCR[0][i][2]) * ipolar + fieldCR.getZ(i);
          } else {
            r[0][i] = 0.0;
            r[1][i] = 0.0;
            r[2][i] = 0.0;
            rCR[0][i] = 0.0;
            rCR[1][i] = 0.0;
            rCR[2][i] = 0.0;
          }

          // Store the current induced dipoles and load the residual induced dipole
          if (use[i]) {
            vec[0][i] = inducedDipole[0][i][0];
            vec[1][i] = inducedDipole[0][i][1];
            vec[2][i] = inducedDipole[0][i][2];
            vecCR[0][i] = inducedDipoleCR[0][i][0];
            vecCR[1][i] = inducedDipoleCR[0][i][1];
            vecCR[2][i] = inducedDipoleCR[0][i][2];
          } else {
            vec[0][i] = 0.0;
            vec[1][i] = 0.0;
            vec[2][i] = 0.0;
            vecCR[0][i] = 0.0;
            vecCR[1][i] = 0.0;
            vecCR[2][i] = 0.0;
          }

          double polar = polarizability[i];
          inducedDipole[0][i][0] = polar * r[0][i];
          inducedDipole[0][i][1] = polar * r[1][i];
          inducedDipole[0][i][2] = polar * r[2][i];
          inducedDipoleCR[0][i][0] = polar * rCR[0][i];
          inducedDipoleCR[0][i][1] = polar * rCR[1][i];
          inducedDipoleCR[0][i][2] = polar * rCR[2][i];
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return IntegerSchedule.fixed();
      }
    }
  }

  /**
   * Set the initial conjugate direction p_0 = z_0 = M^-1 * r_0.
   * <p>
   * Where M^-1 is the inverse of the preconditioner matrix.
   */
  private class InitConjugateRegion extends ParallelRegion {

    private final InitConjugateLoop[] initConjugateLoops;

    public InitConjugateRegion(int nt) {
      initConjugateLoops = new InitConjugateLoop[nt];
    }

    @Override
    public void run() throws Exception {
      try {
        int ti = getThreadIndex();
        if (initConjugateLoops[ti] == null) {
          initConjugateLoops[ti] = new InitConjugateLoop();
        }
        int nAtoms = atoms.length;
        execute(0, nAtoms - 1, initConjugateLoops[ti]);
      } catch (Exception e) {
        String message =
            "Fatal exception computing the mutual induced dipoles in thread "
                + getThreadIndex()
                + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }

    private class InitConjugateLoop extends IntegerForLoop {

      @Override
      public void run(int lb, int ub) throws Exception {

        for (int i = lb; i <= ub; i++) {
          if (use[i]) {
            // Revert to the stored induce dipoles.
            inducedDipole[0][i][0] = vec[0][i];
            inducedDipole[0][i][1] = vec[1][i];
            inducedDipole[0][i][2] = vec[2][i];
            inducedDipoleCR[0][i][0] = vecCR[0][i];
            inducedDipoleCR[0][i][1] = vecCR[1][i];
            inducedDipoleCR[0][i][2] = vecCR[2][i];

            // Set initial conjugate vector p (induced dipoles).
            double polar = polarizability[i];
            z[0][i] = polar * (field.getX(i) + preconditionerScale * r[0][i]);
            z[1][i] = polar * (field.getY(i) + preconditionerScale * r[1][i]);
            z[2][i] = polar * (field.getZ(i) + preconditionerScale * r[2][i]);
            zCR[0][i] = polar * (fieldCR.getX(i) + preconditionerScale * rCR[0][i]);
            zCR[1][i] = polar * (fieldCR.getY(i) + preconditionerScale * rCR[1][i]);
            zCR[2][i] = polar * (fieldCR.getZ(i) + preconditionerScale * rCR[2][i]);
            p[0][i] = z[0][i];
            p[1][i] = z[1][i];
            p[2][i] = z[2][i];
            pCR[0][i] = zCR[0][i];
            pCR[1][i] = zCR[1][i];
            pCR[2][i] = zCR[2][i];
          } else {
            inducedDipole[0][i][0] = 0.0;
            inducedDipole[0][i][1] = 0.0;
            inducedDipole[0][i][2] = 0.0;
            inducedDipoleCR[0][i][0] = 0.0;
            inducedDipoleCR[0][i][1] = 0.0;
            inducedDipoleCR[0][i][2] = 0.0;
            z[0][i] = 0.0;
            z[1][i] = 0.0;
            z[2][i] = 0.0;
            zCR[0][i] = 0.0;
            zCR[1][i] = 0.0;
            zCR[2][i] = 0.0;
            p[0][i] = 0.0;
            p[1][i] = 0.0;
            p[2][i] = 0.0;
            pCR[0][i] = 0.0;
            pCR[1][i] = 0.0;
            pCR[2][i] = 0.0;
          }
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return IntegerSchedule.fixed();
      }
    }
  }

  /**
   * Update the induced dipoles and residuals.
   */
  private class UpdateResidualRegion extends ParallelRegion {

    private final UpdateApLoop[] updateApLoops;
    private final UpdateResidualAndInducedLoop[] updateResidualAndInducedLoops;
    private final SharedDouble pDotApShared;
    private final SharedDouble pDotApCRShared;
    private final SharedDouble rDotZShared;
    private final SharedDouble rDotZCRShared;

    public UpdateResidualRegion(int nt) {
      updateApLoops = new UpdateApLoop[nt];
      updateResidualAndInducedLoops = new UpdateResidualAndInducedLoop[nt];
      pDotApShared = new SharedDouble();
      pDotApCRShared = new SharedDouble();
      rDotZShared = new SharedDouble();
      rDotZCRShared = new SharedDouble();
    }

    public double getRDotZ() {
      return rDotZShared.get();
    }

    public double getRDotZCR() {
      return rDotZCRShared.get();
    }

    @Override
    public void run() throws Exception {
      try {
        int ti = getThreadIndex();
        int nAtoms = atoms.length;
        if (updateApLoops[ti] == null) {
          updateApLoops[ti] = new UpdateApLoop();
          updateResidualAndInducedLoops[ti] = new UpdateResidualAndInducedLoop();
        }
        execute(0, nAtoms - 1, updateApLoops[ti]);
        execute(0, nAtoms - 1, updateResidualAndInducedLoops[ti]);
      } catch (Exception e) {
        String message =
            "Fatal exception computing the mutual induced dipoles in thread "
                + getThreadIndex()
                + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }

    @Override
    public void start() {
      pDotApShared.set(0.0);
      pDotApCRShared.set(0.0);
      rDotZShared.set(0.0);
      rDotZCRShared.set(0.0);
    }

    private class UpdateApLoop extends IntegerForLoop {

      public double pDotAp;
      public double pDotApCR;
      public double rDotZ;
      public double rDotZCR;

      @Override
      public void finish() {
        pDotApShared.addAndGet(pDotAp);
        pDotApCRShared.addAndGet(pDotApCR);
        rDotZShared.addAndGet(rDotZ);
        rDotZCRShared.addAndGet(rDotZCR);
      }

      @Override
      public void run(int lb, int ub) throws Exception {
        for (int i = lb; i <= ub; i++) {
          if (use[i] && polarizability[i] > 0) {
            double ipolar = 1.0 / polarizability[i];
            inducedDipole[0][i][0] = vec[0][i];
            inducedDipole[0][i][1] = vec[1][i];
            inducedDipole[0][i][2] = vec[2][i];
            inducedDipoleCR[0][i][0] = vecCR[0][i];
            inducedDipoleCR[0][i][1] = vecCR[1][i];
            inducedDipoleCR[0][i][2] = vecCR[2][i];
            // Compute Ap = [1/alpha - T] * p
            // Tp is the field due to the conjugate induced dipoles.
            // p/alpha gives the field that induced the conjugate induced dipoles.
            vec[0][i] = p[0][i] * ipolar - field.getX(i);
            vec[1][i] = p[1][i] * ipolar - field.getY(i);
            vec[2][i] = p[2][i] * ipolar - field.getZ(i);
            vecCR[0][i] = pCR[0][i] * ipolar - fieldCR.getX(i);
            vecCR[1][i] = pCR[1][i] * ipolar - fieldCR.getY(i);
            vecCR[2][i] = pCR[2][i] * ipolar - fieldCR.getZ(i);
          } else {
            inducedDipole[0][i][0] = 0.0;
            inducedDipole[0][i][1] = 0.0;
            inducedDipole[0][i][2] = 0.0;
            vec[0][i] = 0.0;
            vec[1][i] = 0.0;
            vec[2][i] = 0.0;
            inducedDipoleCR[0][i][0] = 0.0;
            inducedDipoleCR[0][i][1] = 0.0;
            inducedDipoleCR[0][i][2] = 0.0;
            vecCR[0][i] = 0.0;
            vecCR[1][i] = 0.0;
            vecCR[2][i] = 0.0;
          }

          // Compute dot product of the conjugate vector (p) with Ap (stored in vec).
          pDotAp += p[0][i] * vec[0][i] + p[1][i] * vec[1][i] + p[2][i] * vec[2][i];
          pDotApCR += pCR[0][i] * vecCR[0][i] + pCR[1][i] * vecCR[1][i] + pCR[2][i] * vecCR[2][i];

          // Compute dot product of the residual (r) and preconditioner (z).
          rDotZ += r[0][i] * z[0][i] + r[1][i] * z[1][i] + r[2][i] * z[2][i];
          rDotZCR += rCR[0][i] * zCR[0][i] + rCR[1][i] * zCR[1][i] + rCR[2][i] * zCR[2][i];
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return IntegerSchedule.fixed();
      }

      @Override
      public void start() {
        pDotAp = 0.0;
        pDotApCR = 0.0;
        rDotZ = 0.0;
        rDotZCR = 0.0;
      }
    }

    private class UpdateResidualAndInducedLoop extends IntegerForLoop {

      @Override
      public void run(int lb, int ub) throws Exception {
        double alpha = rDotZShared.get();
        if (pDotApShared.get() != 0) {
          alpha /= pDotApShared.get();
        }
        double alphaCR = rDotZCRShared.get();
        if (pDotApCRShared.get() != 0) {
          alphaCR /= pDotApCRShared.get();
        }
        for (int i = lb; i <= ub; i++) {
          if (use[i]) {
            // Update the residual.
            r[0][i] -= alpha * vec[0][i];
            r[1][i] -= alpha * vec[1][i];
            r[2][i] -= alpha * vec[2][i];
            rCR[0][i] -= alphaCR * vecCR[0][i];
            rCR[1][i] -= alphaCR * vecCR[1][i];
            rCR[2][i] -= alphaCR * vecCR[2][i];

            // Update the induced dipoles based on the scaled conjugate vector.
            // Store them in the work array vec.
            vec[0][i] = inducedDipole[0][i][0] + alpha * p[0][i];
            vec[1][i] = inducedDipole[0][i][1] + alpha * p[1][i];
            vec[2][i] = inducedDipole[0][i][2] + alpha * p[2][i];
            vecCR[0][i] = inducedDipoleCR[0][i][0] + alphaCR * pCR[0][i];
            vecCR[1][i] = inducedDipoleCR[0][i][1] + alphaCR * pCR[1][i];
            vecCR[2][i] = inducedDipoleCR[0][i][2] + alphaCR * pCR[2][i];

            // Prepare to compute the preconditioner.
            // Set the induced dipoles to the polarizability times the residual field.
            double polar = polarizability[i];
            inducedDipole[0][i][0] = polar * r[0][i];
            inducedDipole[0][i][1] = polar * r[1][i];
            inducedDipole[0][i][2] = polar * r[2][i];
            inducedDipoleCR[0][i][0] = polar * rCR[0][i];
            inducedDipoleCR[0][i][1] = polar * rCR[1][i];
            inducedDipoleCR[0][i][2] = polar * rCR[2][i];
          } else {
            vec[0][i] = 0.0;
            vec[1][i] = 0.0;
            vec[2][i] = 0.0;
            vecCR[0][i] = 0.0;
            vecCR[1][i] = 0.0;
            vecCR[2][i] = 0.0;
            inducedDipole[0][i][0] = 0.0;
            inducedDipole[0][i][1] = 0.0;
            inducedDipole[0][i][2] = 0.0;
            inducedDipoleCR[0][i][0] = 0.0;
            inducedDipoleCR[0][i][1] = 0.0;
            inducedDipoleCR[0][i][2] = 0.0;
          }
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return IntegerSchedule.fixed();
      }
    }
  }

  /**
   * Update the preconditioner and conjugate search direction.
   */
  private class UpdateConjugateRegion extends ParallelRegion {

    private final UpdatePreconditionerLoop[] updatePreconditionerLoops;
    private final UpdateConjugateLoop[] updateConjugateLoops;
    private final SharedDouble betaShared;
    private final SharedDouble betaCRShared;
    private final SharedDouble epsShared;
    private final SharedDouble epsCRShared;
    public double previousRDotZ;
    public double previousRDotZCR;

    public UpdateConjugateRegion(int nt) {
      updatePreconditionerLoops = new UpdatePreconditionerLoop[nt];
      updateConjugateLoops = new UpdateConjugateLoop[nt];
      betaShared = new SharedDouble();
      betaCRShared = new SharedDouble();
      epsShared = new SharedDouble();
      epsCRShared = new SharedDouble();
    }

    public double getEps() {
      return epsShared.get();
    }

    public double getEpsCR() {
      return epsCRShared.get();
    }

    @Override
    public void run() throws Exception {
      try {
        int ti = getThreadIndex();
        if (updatePreconditionerLoops[ti] == null) {
          updatePreconditionerLoops[ti] = new UpdatePreconditionerLoop();
          updateConjugateLoops[ti] = new UpdateConjugateLoop();
        }
        int nAtoms = atoms.length;
        execute(0, nAtoms - 1, updatePreconditionerLoops[ti]);
        execute(0, nAtoms - 1, updateConjugateLoops[ti]);
      } catch (Exception e) {
        String message =
            "Fatal exception computing the mutual induced dipoles in thread "
                + getThreadIndex()
                + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }

    @Override
    public void start() {
      betaShared.set(0.0);
      betaCRShared.set(0.0);
      epsShared.set(0.0);
      epsCRShared.set(0.0);
      if (previousRDotZ == 0.0) {
        previousRDotZ = 1.0;
      }
      if (previousRDotZCR == 0.0) {
        previousRDotZCR = 1.0;
      }
    }

    private class UpdatePreconditionerLoop extends IntegerForLoop {

      public double rDotZ;
      public double rDotZCR;

      @Override
      public void finish() {
        betaShared.addAndGet(rDotZ / previousRDotZ);
        betaCRShared.addAndGet(rDotZCR / previousRDotZCR);
      }

      @Override
      public void run(int lb, int ub) throws Exception {
        for (int i = lb; i <= ub; i++) {
          if (use[i]) {
            // Revert the induced dipoles to the saved values.
            inducedDipole[0][i][0] = vec[0][i];
            inducedDipole[0][i][1] = vec[1][i];
            inducedDipole[0][i][2] = vec[2][i];
            inducedDipoleCR[0][i][0] = vecCR[0][i];
            inducedDipoleCR[0][i][1] = vecCR[1][i];
            inducedDipoleCR[0][i][2] = vecCR[2][i];

            // Compute the dot product of the residual and preconditioner.
            double polar = polarizability[i];
            z[0][i] = polar * (field.getX(i) + preconditionerScale * r[0][i]);
            z[1][i] = polar * (field.getY(i) + preconditionerScale * r[1][i]);
            z[2][i] = polar * (field.getZ(i) + preconditionerScale * r[2][i]);
            zCR[0][i] = polar * (fieldCR.getX(i) + preconditionerScale * rCR[0][i]);
            zCR[1][i] = polar * (fieldCR.getY(i) + preconditionerScale * rCR[1][i]);
            zCR[2][i] = polar * (fieldCR.getZ(i) + preconditionerScale * rCR[2][i]);
            rDotZ += r[0][i] * z[0][i] + r[1][i] * z[1][i] + r[2][i] * z[2][i];
            rDotZCR += rCR[0][i] * zCR[0][i] + rCR[1][i] * zCR[1][i] + rCR[2][i] * zCR[2][i];
          } else {
            inducedDipole[0][i][0] = 0.0;
            inducedDipole[0][i][1] = 0.0;
            inducedDipole[0][i][2] = 0.0;
            inducedDipoleCR[0][i][0] = 0.0;
            inducedDipoleCR[0][i][1] = 0.0;
            inducedDipoleCR[0][i][2] = 0.0;
            z[0][i] = 0.0;
            z[1][i] = 0.0;
            z[2][i] = 0.0;
            zCR[0][i] = 0.0;
            zCR[1][i] = 0.0;
            zCR[2][i] = 0.0;
          }
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return IntegerSchedule.fixed();
      }

      @Override
      public void start() {
        rDotZ = 0.0;
        rDotZCR = 0.0;
      }
    }

    private class UpdateConjugateLoop extends IntegerForLoop {

      public double eps;
      public double epsCR;

      @Override
      public void finish() {
        epsShared.addAndGet(eps);
        epsCRShared.addAndGet(epsCR);
      }

      @Override
      public void run(int lb, int ub) throws Exception {
        double beta = betaShared.get();
        double betaCR = betaCRShared.get();
        for (int i = lb; i <= ub; i++) {
          if (use[i]) {
            // Update the conjugate vector.
            p[0][i] = z[0][i] + beta * p[0][i];
            p[1][i] = z[1][i] + beta * p[1][i];
            p[2][i] = z[2][i] + beta * p[2][i];
            pCR[0][i] = zCR[0][i] + betaCR * pCR[0][i];
            pCR[1][i] = zCR[1][i] + betaCR * pCR[1][i];
            pCR[2][i] = zCR[2][i] + betaCR * pCR[2][i];
            // Sum the square of the residual field.
            eps += r[0][i] * r[0][i] + r[1][i] * r[1][i] + r[2][i] * r[2][i];
            epsCR += rCR[0][i] * rCR[0][i] + rCR[1][i] * rCR[1][i] + rCR[2][i] * rCR[2][i];
          } else {
            p[0][i] = 0.0;
            p[1][i] = 0.0;
            p[2][i] = 0.0;
            pCR[0][i] = 0.0;
            pCR[1][i] = 0.0;
            pCR[2][i] = 0.0;
          }
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return IntegerSchedule.fixed();
      }

      @Override
      public void start() {
        eps = 0.0;
        epsCR = 0.0;
      }
    }
  }

  /** Evaluate the real space field due to induced dipoles using a short cutoff (~3-4 A). */
  private class PreconditionerRegion extends ParallelRegion {

    private final InducedPreconditionerFieldLoop[] inducedPreconditionerFieldLoop;

    PreconditionerRegion(int threadCount) {
      inducedPreconditionerFieldLoop = new InducedPreconditionerFieldLoop[threadCount];
    }

    @Override
    public void run() {
      int threadIndex = getThreadIndex();
      if (inducedPreconditionerFieldLoop[threadIndex] == null) {
        inducedPreconditionerFieldLoop[threadIndex] = new InducedPreconditionerFieldLoop();
      }
      try {
        int nAtoms = atoms.length;
        execute(0, nAtoms - 1, inducedPreconditionerFieldLoop[threadIndex]);
      } catch (Exception e) {
        String message =
            "Fatal exception computing the induced real space field in thread "
                + getThreadIndex()
                + "\n";
        logger.log(Level.SEVERE, message, e);
      }
    }

    private class InducedPreconditionerFieldLoop extends IntegerForLoop {

      private int threadID;
      private double[] x, y, z;
      private double[][] ind, indCR;

      InducedPreconditionerFieldLoop() {
      }

      @Override
      public void finish() {
        realSpaceSCFTime[threadID] += System.nanoTime();
      }

      @Override
      public void run(int lb, int ub) {
        final double[] dx = new double[3];
        final double[][] transOp = new double[3][3];

        // Loop over a chunk of atoms.
        int[][] lists = preconditionerLists[0];
        int[] counts = preconditionerCounts[0];
        for (int i = lb; i <= ub; i++) {
          if (!use[i]) {
            continue;
          }
          double fx = 0.0;
          double fy = 0.0;
          double fz = 0.0;
          double px = 0.0;
          double py = 0.0;
          double pz = 0.0;
          final double xi = x[i];
          final double yi = y[i];
          final double zi = z[i];
          final double[] dipolei = ind[i];
          final double uix = dipolei[0];
          final double uiy = dipolei[1];
          final double uiz = dipolei[2];
          final double[] dipoleCRi = indCR[i];
          final double pix = dipoleCRi[0];
          final double piy = dipoleCRi[1];
          final double piz = dipoleCRi[2];
          final double pdi = ipdamp[i];
          final double pti = thole[i];

          // Loop over the neighbor list.
          final int[] list = lists[i];
          final int npair = counts[i];
          for (int j = 0; j < npair; j++) {
            final int k = list[j];
            if (!use[k]) {
              continue;
            }
            final double pdk = ipdamp[k];
            final double ptk = thole[k];
            dx[0] = x[k] - xi;
            dx[1] = y[k] - yi;
            dx[2] = z[k] - zi;
            final double r2 = crystal.image(dx);

            // Calculate the error function damping terms.
            final double r = sqrt(r2);
            final double rr1 = 1.0 / r;
            final double rr2 = rr1 * rr1;
            final double ralpha = ewaldParameters.aewald * r;
            final double exp2a = exp(-ralpha * ralpha);
            final double bn0 = erfc(ralpha) * rr1;
            // final double exp2a = 1.0;
            // final double bn0 = rr1;
            final double bn1 = (bn0 + ewaldParameters.an0 * exp2a) * rr2;
            final double bn2 = (3.0 * bn1 + ewaldParameters.an1 * exp2a) * rr2;
            double scale3 = 1.0;
            double scale5 = 1.0;
            double damp = pdi * pdk;
            final double pgamma = min(pti, ptk);
            final double rdamp = r * damp;
            damp = -pgamma * rdamp * rdamp * rdamp;
            if (damp > -50.0) {
              final double expdamp = exp(damp);
              scale3 = 1.0 - expdamp;
              scale5 = 1.0 - expdamp * (1.0 - damp);
            }
            double rr3 = rr1 * rr2;
            double rr5 = 3.0 * rr3 * rr2;
            rr3 *= (1.0 - scale3);
            rr5 *= (1.0 - scale5);
            final double xr = dx[0];
            final double yr = dx[1];
            final double zr = dx[2];
            final double[] dipolek = ind[k];
            final double ukx = dipolek[0];
            final double uky = dipolek[1];
            final double ukz = dipolek[2];
            final double ukr = ukx * xr + uky * yr + ukz * zr;
            final double bn2ukr = bn2 * ukr;
            final double fimx = -bn1 * ukx + bn2ukr * xr;
            final double fimy = -bn1 * uky + bn2ukr * yr;
            final double fimz = -bn1 * ukz + bn2ukr * zr;
            final double rr5ukr = rr5 * ukr;
            final double fidx = -rr3 * ukx + rr5ukr * xr;
            final double fidy = -rr3 * uky + rr5ukr * yr;
            final double fidz = -rr3 * ukz + rr5ukr * zr;
            fx += (fimx - fidx);
            fy += (fimy - fidy);
            fz += (fimz - fidz);
            final double[] dipolepk = indCR[k];
            final double pkx = dipolepk[0];
            final double pky = dipolepk[1];
            final double pkz = dipolepk[2];
            final double pkr = pkx * xr + pky * yr + pkz * zr;
            final double bn2pkr = bn2 * pkr;
            final double pimx = -bn1 * pkx + bn2pkr * xr;
            final double pimy = -bn1 * pky + bn2pkr * yr;
            final double pimz = -bn1 * pkz + bn2pkr * zr;
            final double rr5pkr = rr5 * pkr;
            final double pidx = -rr3 * pkx + rr5pkr * xr;
            final double pidy = -rr3 * pky + rr5pkr * yr;
            final double pidz = -rr3 * pkz + rr5pkr * zr;
            px += (pimx - pidx);
            py += (pimy - pidy);
            pz += (pimz - pidz);
            final double uir = uix * xr + uiy * yr + uiz * zr;
            final double bn2uir = bn2 * uir;
            final double fkmx = -bn1 * uix + bn2uir * xr;
            final double fkmy = -bn1 * uiy + bn2uir * yr;
            final double fkmz = -bn1 * uiz + bn2uir * zr;
            final double rr5uir = rr5 * uir;
            final double fkdx = -rr3 * uix + rr5uir * xr;
            final double fkdy = -rr3 * uiy + rr5uir * yr;
            final double fkdz = -rr3 * uiz + rr5uir * zr;
            field.add(threadID, k, fkmx - fkdx, fkmy - fkdy, fkmz - fkdz);
            final double pir = pix * xr + piy * yr + piz * zr;
            final double bn2pir = bn2 * pir;
            final double pkmx = -bn1 * pix + bn2pir * xr;
            final double pkmy = -bn1 * piy + bn2pir * yr;
            final double pkmz = -bn1 * piz + bn2pir * zr;
            final double rr5pir = rr5 * pir;
            final double pkdx = -rr3 * pix + rr5pir * xr;
            final double pkdy = -rr3 * piy + rr5pir * yr;
            final double pkdz = -rr3 * piz + rr5pir * zr;
            fieldCR.add(threadID, k, pkmx - pkdx, pkmy - pkdy, pkmz - pkdz);
          }
          field.add(threadID, i, fx, fy, fz);
          fieldCR.add(threadID, i, px, py, pz);
        }

        // Loop over symmetry mates.
        List<SymOp> symOps = crystal.spaceGroup.symOps;
        int nSymm = symOps.size();
        for (int iSymm = 1; iSymm < nSymm; iSymm++) {
          SymOp symOp = crystal.spaceGroup.getSymOp(iSymm);
          crystal.getTransformationOperator(symOp, transOp);
          lists = preconditionerLists[iSymm];
          counts = preconditionerCounts[iSymm];
          final double[] xs = coordinates[iSymm][0];
          final double[] ys = coordinates[iSymm][1];
          final double[] zs = coordinates[iSymm][2];
          final double[][] inds = inducedDipole[iSymm];
          final double[][] indCRs = inducedDipoleCR[iSymm];

          // Loop over a chunk of atoms.
          for (int i = lb; i <= ub; i++) {
            if (!use[i]) {
              continue;
            }
            double fx = 0.0;
            double fy = 0.0;
            double fz = 0.0;
            double px = 0.0;
            double py = 0.0;
            double pz = 0.0;
            final double xi = x[i];
            final double yi = y[i];
            final double zi = z[i];
            final double[] dipolei = ind[i];
            final double uix = dipolei[0];
            final double uiy = dipolei[1];
            final double uiz = dipolei[2];
            final double[] dipoleCRi = indCR[i];
            final double pix = dipoleCRi[0];
            final double piy = dipoleCRi[1];
            final double piz = dipoleCRi[2];
            final double pdi = ipdamp[i];
            final double pti = thole[i];

            // Loop over the neighbor list.
            final int[] list = lists[i];
            final int npair = counts[i];
            for (int j = 0; j < npair; j++) {
              final int k = list[j];
              if (!use[k]) {
                continue;
              }
              double selfScale = 1.0;
              if (i == k) {
                selfScale = 0.5;
              }
              final double pdk = ipdamp[k];
              final double ptk = thole[k];
              dx[0] = xs[k] - xi;
              dx[1] = ys[k] - yi;
              dx[2] = zs[k] - zi;
              final double r2 = crystal.image(dx);

              // Calculate the error function damping terms.
              final double r = sqrt(r2);
              final double rr1 = 1.0 / r;
              final double rr2 = rr1 * rr1;
              final double ralpha = ewaldParameters.aewald * r;
              final double exp2a = exp(-ralpha * ralpha);
              final double bn0 = erfc(ralpha) * rr1;
              // final double exp2a = 1.0;
              // final double bn0 = rr1;
              final double bn1 = (bn0 + ewaldParameters.an0 * exp2a) * rr2;
              final double bn2 = (3.0 * bn1 + ewaldParameters.an1 * exp2a) * rr2;
              double scale3 = 1.0;
              double scale5 = 1.0;
              double damp = pdi * pdk;
              final double pgamma = min(pti, ptk);
              final double rdamp = r * damp;
              damp = -pgamma * rdamp * rdamp * rdamp;
              if (damp > -50.0) {
                final double expdamp = exp(damp);
                scale3 = 1.0 - expdamp;
                scale5 = 1.0 - expdamp * (1.0 - damp);
              }
              double rr3 = rr1 * rr2;
              double rr5 = 3.0 * rr3 * rr2;
              rr3 *= (1.0 - scale3);
              rr5 *= (1.0 - scale5);
              final double xr = dx[0];
              final double yr = dx[1];
              final double zr = dx[2];
              final double[] dipolek = inds[k];
              final double ukx = dipolek[0];
              final double uky = dipolek[1];
              final double ukz = dipolek[2];
              final double[] dipolepk = indCRs[k];
              final double pkx = dipolepk[0];
              final double pky = dipolepk[1];
              final double pkz = dipolepk[2];
              final double ukr = ukx * xr + uky * yr + ukz * zr;
              final double bn2ukr = bn2 * ukr;
              final double fimx = -bn1 * ukx + bn2ukr * xr;
              final double fimy = -bn1 * uky + bn2ukr * yr;
              final double fimz = -bn1 * ukz + bn2ukr * zr;
              final double rr5ukr = rr5 * ukr;
              final double fidx = -rr3 * ukx + rr5ukr * xr;
              final double fidy = -rr3 * uky + rr5ukr * yr;
              final double fidz = -rr3 * ukz + rr5ukr * zr;
              fx += selfScale * (fimx - fidx);
              fy += selfScale * (fimy - fidy);
              fz += selfScale * (fimz - fidz);
              final double pkr = pkx * xr + pky * yr + pkz * zr;
              final double bn2pkr = bn2 * pkr;
              final double pimx = -bn1 * pkx + bn2pkr * xr;
              final double pimy = -bn1 * pky + bn2pkr * yr;
              final double pimz = -bn1 * pkz + bn2pkr * zr;
              final double rr5pkr = rr5 * pkr;
              final double pidx = -rr3 * pkx + rr5pkr * xr;
              final double pidy = -rr3 * pky + rr5pkr * yr;
              final double pidz = -rr3 * pkz + rr5pkr * zr;
              px += selfScale * (pimx - pidx);
              py += selfScale * (pimy - pidy);
              pz += selfScale * (pimz - pidz);
              final double uir = uix * xr + uiy * yr + uiz * zr;
              final double bn2uir = bn2 * uir;
              final double fkmx = -bn1 * uix + bn2uir * xr;
              final double fkmy = -bn1 * uiy + bn2uir * yr;
              final double fkmz = -bn1 * uiz + bn2uir * zr;
              final double rr5uir = rr5 * uir;
              final double fkdx = -rr3 * uix + rr5uir * xr;
              final double fkdy = -rr3 * uiy + rr5uir * yr;
              final double fkdz = -rr3 * uiz + rr5uir * zr;
              double xc = selfScale * (fkmx - fkdx);
              double yc = selfScale * (fkmy - fkdy);
              double zc = selfScale * (fkmz - fkdz);
              double fkx = (xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0]);
              double fky = (xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1]);
              double fkz = (xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2]);
              field.add(threadID, k, fkx, fky, fkz);
              final double pir = pix * xr + piy * yr + piz * zr;
              final double bn2pir = bn2 * pir;
              final double pkmx = -bn1 * pix + bn2pir * xr;
              final double pkmy = -bn1 * piy + bn2pir * yr;
              final double pkmz = -bn1 * piz + bn2pir * zr;
              final double rr5pir = rr5 * pir;
              final double pkdx = -rr3 * pix + rr5pir * xr;
              final double pkdy = -rr3 * piy + rr5pir * yr;
              final double pkdz = -rr3 * piz + rr5pir * zr;
              xc = selfScale * (pkmx - pkdx);
              yc = selfScale * (pkmy - pkdy);
              zc = selfScale * (pkmz - pkdz);
              fkx = (xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0]);
              fky = (xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1]);
              fkz = (xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2]);
              fieldCR.add(threadID, k, fkx, fky, fkz);
            }
            field.add(threadID, i, fx, fy, fz);
            fieldCR.add(threadID, i, px, py, pz);
          }
        }
      }

      @Override
      public IntegerSchedule schedule() {
        return realSpaceSchedule;
      }

      @Override
      public void start() {
        threadID = getThreadIndex();
        realSpaceSCFTime[threadID] -= System.nanoTime();
        x = coordinates[0][0];
        y = coordinates[0][1];
        z = coordinates[0][2];
        ind = inducedDipole[0];
        indCR = inducedDipoleCR[0];
      }
    }
  }
}
