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

import static ffx.potential.parameters.MultipoleType.getRotationMatrix;
import static ffx.potential.parameters.MultipoleType.rotateMultipole;
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
import static org.apache.commons.math3.util.FastMath.max;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import ffx.potential.parameters.PolarizeType;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Parallel initialization of accumulation arrays, expand atomic coordinates and rotation of
 * multipoles into the global frame.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class InitializationRegion extends ParallelRegion {

  private static final Logger logger = Logger.getLogger(InitializationRegion.class.getName());

  private final ParticleMeshEwald particleMeshEwald;
  /**
   * If set to false, multipoles are fixed in their local frame and torques are zero, which is useful
   * for narrowing down discrepancies between analytic and finite-difference derivatives(default is
   * true).
   */
  private final boolean rotateMultipoles;
  /** If set to false, multipole charges are set to zero (default is true). */
  private final boolean useCharges;
  /** If set to false, multipole dipoles are set to zero (default is true). */
  private final boolean useDipoles;
  /** If set to false, multipole quadrupoles are set to zero (default is true). */
  private final boolean useQuadrupoles;
  /**
   * Initialization Loops
   */
  private final InitializationLoop[] initializationLoop;
  /**
   * Rotate Multipole Loops
   */
  private final RotateMultipolesLoop[] rotateMultipolesLoop;
  /**
   * If lambdaTerm is true, some ligand atom interactions with the environment are being turned
   * on/off.
   */
  private boolean lambdaTerm;

  private ExtendedSystem esvSystem;
  /**
   * If esvTerm is true, the electrostatics of some atoms are being titrated.
   */
  private boolean esvTerm;

  /** An ordered array of atoms in the system. */
  private Atom[] atoms;
  /** Dimensions of [nsymm][xyz][nAtoms]. */
  private double[][][] coordinates;
  /** Unit cell and spacegroup information. */
  private Crystal crystal;
  /** Multipole frame definition. */
  private MultipoleFrameDefinition[] frame;
  /** Multipole frame defining atoms. */
  private int[][] axisAtom;
  /** Dimensions of [nsymm][nAtoms][10] */
  private double[][][] globalMultipole;
  /** Dimensions of [nsymm][nAtoms][10] */
  private double[][][] titrationMultipole;
  private double[][][] tautomerMultipole;
  /** Polarizability of each atom */
  private double[] polarizability;
  private double[] titrationPolarizability;
  private double[] tautomerPolarizability;
  private double[] thole;
  private double[] ipdamp;
  /**
   * The "use" array can be employed to turn off atoms for computing the electrostatic energy of
   * sub-structures.
   */
  private boolean[] use;
  /**
   * Neighbor lists, including atoms beyond the real space cutoff. [nsymm][nAtoms][nAllNeighbors]
   */
  private int[][][] neighborLists;
  /**
   * Neighbor lists, without atoms beyond the real space cutoff. [nSymm][nAtoms][nIncludedNeighbors]
   */
  private int[][][] realSpaceLists;

  private int[][][] vaporLists;
  /** Atomic Gradient array. */
  private AtomicDoubleArray3D grad;
  /** Atomic Torque array. */
  private AtomicDoubleArray3D torque;
  /** Partial derivative of the gradient with respect to Lambda. */
  private AtomicDoubleArray3D lambdaGrad;
  /** Partial derivative of the torque with respect to Lambda. */
  private AtomicDoubleArray3D lambdaTorque;

  public InitializationRegion(ParticleMeshEwald particleMeshEwald, int maxThreads,
      ForceField forceField) {
    initializationLoop = new InitializationLoop[maxThreads];
    rotateMultipolesLoop = new RotateMultipolesLoop[maxThreads];
    useCharges = forceField.getBoolean("USE_CHARGES", true);
    useDipoles = forceField.getBoolean("USE_DIPOLES", true);
    useQuadrupoles = forceField.getBoolean("USE_QUADRUPOLES", true);
    rotateMultipoles = forceField.getBoolean("ROTATE_MULTIPOLES", true);
    this.particleMeshEwald = particleMeshEwald;
  }

  /**
   * Execute the InitializationRegion with the passed ParallelTeam.
   *
   * @param parallelTeam The ParallelTeam instance to execute with.
   */
  public void executeWith(ParallelTeam parallelTeam) {
    try {
      parallelTeam.execute(this);
    } catch (RuntimeException e) {
      String message = "RuntimeException expanding coordinates and rotating multipoles.\n";
      logger.log(Level.WARNING, message, e);
      throw e;
    } catch (Exception e) {
      String message = "Fatal exception expanding coordinates and rotating multipoles.\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  public void init(
      boolean lambdaTerm,
      ExtendedSystem esvSystem,
      Atom[] atoms,
      double[][][] coordinates,
      Crystal crystal,
      MultipoleFrameDefinition[] frame,
      int[][] axisAtom,
      double[][][] globalMultipole,
      double[][][] titrationMultipole,
      double[][][] tautomerMultipole,
      double[] polarizability,
      double[] titrationPolarizability,
      double[] tautomerPolarizability,
      double[] thole,
      double[] ipdamp,
      boolean[] use,
      int[][][] neighborLists,
      int[][][] realSpaceLists,
      int[][][] vaporLists,
      AtomicDoubleArray3D grad,
      AtomicDoubleArray3D torque,
      AtomicDoubleArray3D lambdaGrad,
      AtomicDoubleArray3D lambdaTorque) {
    this.lambdaTerm = lambdaTerm;
    this.esvSystem = esvSystem;
    if (esvSystem != null) {
      this.esvTerm = true;
    }
    this.atoms = atoms;
    this.coordinates = coordinates;
    this.crystal = crystal;
    this.frame = frame;
    this.axisAtom = axisAtom;
    this.globalMultipole = globalMultipole;
    this.titrationMultipole = titrationMultipole;
    this.tautomerMultipole = tautomerMultipole;
    this.polarizability = polarizability;
    this.titrationPolarizability = titrationPolarizability;
    this.tautomerPolarizability = tautomerPolarizability;
    this.thole = thole;
    this.ipdamp = ipdamp;
    this.use = use;
    this.neighborLists = neighborLists;
    this.realSpaceLists = realSpaceLists;
    this.vaporLists = vaporLists;
    this.grad = grad;
    this.torque = torque;
    this.lambdaGrad = lambdaGrad;
    this.lambdaTorque = lambdaTorque;
  }

  @Override
  public void run() {

    int nAtoms = atoms.length;
    int threadIndex = getThreadIndex();
    if (initializationLoop[threadIndex] == null) {
      initializationLoop[threadIndex] = new InitializationLoop();
      rotateMultipolesLoop[threadIndex] = new RotateMultipolesLoop();
    }
    try {
      execute(0, nAtoms - 1, initializationLoop[threadIndex]);
      execute(0, nAtoms - 1, rotateMultipolesLoop[threadIndex]);
    } catch (Exception e) {
      String message = "Fatal exception initializing coordinates in thread: " + threadIndex + "\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  private class InitializationLoop extends IntegerForLoop {

    private final double[] in = new double[3];
    private final double[] out = new double[3];
    private double[] x;
    private double[] y;
    private double[] z;
    private int threadID;
    // Extra padding to avert cache interference.
    private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
    private long pad8, pad9, pada, padb, padc, padd, pade, padf;

    @Override
    public void run(int lb, int ub) {
      grad.reset(threadID, lb, ub);
      torque.reset(threadID, lb, ub);
      if (lambdaTerm) {
        lambdaGrad.reset(threadID, lb, ub);
        lambdaTorque.reset(threadID, lb, ub);
      }

      // Initialize the local coordinate arrays.
      for (int i = lb; i <= ub; i++) {
        Atom atom = atoms[i];
        x[i] = atom.getX();
        y[i] = atom.getY();
        z[i] = atom.getZ();
        use[i] = atom.getUse();
        /*
         Real space Ewald is cutoff at ~7 A, compared to ~12 A for vdW,
         so the number of neighbors is much more compact. A specific list for real space Ewald is filled during
         computation of the permanent real space field that includes only evaluated interactions. Subsequent real
         space loops, especially the SCF, then do not spend time evaluating pairwise distances outside the cutoff.
        */
        int size = neighborLists[0][i].length;
        if (vaporLists != null) {
          size = max(size, vaporLists[0][i].length);
        }
        if (realSpaceLists[0][i] == null || realSpaceLists[0][i].length < size) {
          realSpaceLists[0][i] = new int[size];
        }
      }

      // Expand coordinates.
      List<SymOp> symOps = crystal.spaceGroup.symOps;
      int nSymm = symOps.size();
      for (int iSymm = 1; iSymm < nSymm; iSymm++) {
        SymOp symOp = symOps.get(iSymm);
        double[] xs = coordinates[iSymm][0];
        double[] ys = coordinates[iSymm][1];
        double[] zs = coordinates[iSymm][2];
        for (int i = lb; i <= ub; i++) {
          in[0] = x[i];
          in[1] = y[i];
          in[2] = z[i];
          crystal.applySymOp(in, out, symOp);
          xs[i] = out[0];
          ys[i] = out[1];
          zs[i] = out[2];
          int size = neighborLists[iSymm][i].length;
          if (realSpaceLists[iSymm][i] == null || realSpaceLists[iSymm][i].length < size) {
            realSpaceLists[iSymm][i] = new int[size];
          }
        }
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return IntegerSchedule.fixed();
    }

    @Override
    public void start() {
      x = coordinates[0][0];
      y = coordinates[0][1];
      z = coordinates[0][2];
      threadID = getThreadIndex();
    }
  }

  private class RotateMultipolesLoop extends IntegerForLoop {

    // Local variables
    private final double[] localOrigin = new double[3];
    private final double[][] frameCoords = new double[4][3];
    private final double[][] rotmat = new double[3][3];
    private final double[] tempDipole = new double[3];
    private final double[][] tempQuadrupole = new double[3][3];
    private final double[] dipole = new double[3];
    private final double[][] quadrupole = new double[3][3];
    // Extra padding to avert cache interference.
    private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
    private long pad8, pad9, pada, padb, padc, padd, pade, padf;

    @Override
    public void run(int lb, int ub) {
      List<SymOp> symOps = crystal.spaceGroup.symOps;
      int nSymm = symOps.size();
      for (int iSymm = 0; iSymm < nSymm; iSymm++) {
        final double[] x = coordinates[iSymm][0];
        final double[] y = coordinates[iSymm][1];
        final double[] z = coordinates[iSymm][2];
        for (int ii = lb; ii <= ub; ii++) {
          Atom atom = atoms[ii];
          double chargeScale = 1.0;
          double dipoleScale = 1.0;
          double quadrupoleScale = 1.0;
          double polarizabilityScale = 1.0;
          if (!useCharges) {
            chargeScale = 0.0;
          }
          if (!useDipoles) {
            dipoleScale = 0.0;
          }
          if (!useQuadrupoles) {
            quadrupoleScale = 0.0;
          }

          double elecScale = 1.0;
          if (!atom.getElectrostatics()) {
            elecScale = 0.0;
          }

          // Collect the MultipoleType for Atom i.
          MultipoleType multipoleType = particleMeshEwald.getMultipoleType(ii);
          double[] in = multipoleType.getMultipole();
          // Update the frame
          frame[ii] = multipoleType.frameDefinition;
          // Update the axis defining atom.
          axisAtom[ii] = atom.getAxisAtomIndices();

          if (rotateMultipoles) {
            // Local frame origin is the location of the current atomic multipole atom.
            localOrigin[0] = x[ii];
            localOrigin[1] = y[ii];
            localOrigin[2] = z[ii];

            // Collect coordinates of the frame defining atoms.
            int[] referenceSites = axisAtom[ii];
            int nSites = 0;
            if (referenceSites != null) {
              nSites = referenceSites.length;
            }
            for (int i = 0; i < nSites; i++) {
              int index = referenceSites[i];
              frameCoords[i][0] = x[index];
              frameCoords[i][1] = y[index];
              frameCoords[i][2] = z[index];
            }

            // Load the dipole for rotation.
            tempDipole[0] = in[t100];
            tempDipole[1] = in[t010];
            tempDipole[2] = in[t001];
            // Load the quadrupole for rotation.
            tempQuadrupole[0][0] = in[t200];
            tempQuadrupole[1][1] = in[t020];
            tempQuadrupole[2][2] = in[t002];
            tempQuadrupole[0][1] = in[t110];
            tempQuadrupole[0][2] = in[t101];
            tempQuadrupole[1][2] = in[t011];
            tempQuadrupole[1][0] = in[t110];
            tempQuadrupole[2][0] = in[t101];
            tempQuadrupole[2][1] = in[t011];
            // Check for chiral flipping.

            boolean needsChiralInversion = false;
            /*
            boolean needsChiralInversion = checkMultipoleChirality(frame[ii], localOrigin, frameCoords);
            if (needsChiralInversion) {
              // Flip the sign of the Y-dipole
              tempDipole[1] = -tempDipole[1];
              // Flip the sign of the XY-quadrupole
              tempQuadrupole[0][1] = -tempQuadrupole[0][1];
              tempQuadrupole[1][0] = -tempQuadrupole[1][0];
              // Flip the sign of the YZ-quadrupole
              tempQuadrupole[1][2] = -tempQuadrupole[1][2];
              tempQuadrupole[2][1] = -tempQuadrupole[2][1];
            }
            */
            getRotationMatrix(frame[ii], localOrigin, frameCoords, rotmat);
            rotateMultipole(rotmat, tempDipole, tempQuadrupole, dipole, quadrupole);

            double[] out = globalMultipole[iSymm][ii];
            // Set the charge.
            out[t000] = in[0] * chargeScale * elecScale;
            // Set the dipole in the global frame.
            out[t100] = dipole[0] * dipoleScale * elecScale;
            out[t010] = dipole[1] * dipoleScale * elecScale;
            out[t001] = dipole[2] * dipoleScale * elecScale;
            // Set the quadrupole in the global frame.
            out[t200] = quadrupole[0][0] * quadrupoleScale * elecScale;
            out[t020] = quadrupole[1][1] * quadrupoleScale * elecScale;
            out[t002] = quadrupole[2][2] * quadrupoleScale * elecScale;
            out[t110] = quadrupole[0][1] * quadrupoleScale * elecScale;
            out[t101] = quadrupole[0][2] * quadrupoleScale * elecScale;
            out[t011] = quadrupole[1][2] * quadrupoleScale * elecScale;
            /* For ESV atoms, also rotate and scale the Mdot multipole. */
            if (esvTerm && esvSystem.isTitrating(ii)) {
              double[] esvMultipoleTitrDot = new double[10];
              double[] esvMultipoleTautDot = new double[10];
              double titrationLambda = esvSystem.getTitrationLambda(ii);
              double tautomerLambda = esvSystem.getTautomerLambda(ii);
              esvMultipoleTitrDot = esvSystem.getTitrationUtils()
                  .getMultipoleTitrationDeriv(atom, titrationLambda,
                      tautomerLambda, esvMultipoleTitrDot);
              esvMultipoleTautDot = esvSystem.getTitrationUtils()
                  .getMultipoleTautomerDeriv(atom, titrationLambda,
                      tautomerLambda, esvMultipoleTautDot);
              double tempCharge = esvMultipoleTitrDot[0];
              // Load the dipole for rotation.
              tempDipole[0] = esvMultipoleTitrDot[t100];
              tempDipole[1] = esvMultipoleTitrDot[t010];
              tempDipole[2] = esvMultipoleTitrDot[t001];
              // Load the quadrupole for rotation.
              tempQuadrupole[0][0] = esvMultipoleTitrDot[t200];
              tempQuadrupole[1][1] = esvMultipoleTitrDot[t020];
              tempQuadrupole[2][2] = esvMultipoleTitrDot[t002];
              tempQuadrupole[0][1] = esvMultipoleTitrDot[t110];
              tempQuadrupole[0][2] = esvMultipoleTitrDot[t101];
              tempQuadrupole[1][2] = esvMultipoleTitrDot[t011];
              tempQuadrupole[1][0] = esvMultipoleTitrDot[t110];
              tempQuadrupole[2][0] = esvMultipoleTitrDot[t101];
              tempQuadrupole[2][1] = esvMultipoleTitrDot[t011];
              if (needsChiralInversion) {
                tempDipole[1] = -tempDipole[1];
                tempQuadrupole[0][1] = -tempQuadrupole[0][1];
                tempQuadrupole[1][0] = -tempQuadrupole[1][0];
                tempQuadrupole[1][2] = -tempQuadrupole[1][2];
                tempQuadrupole[2][1] = -tempQuadrupole[2][1];
              }
              rotateMultipole(rotmat, tempDipole, tempQuadrupole, dipole, quadrupole);
              out = titrationMultipole[iSymm][ii];
              out[t000] = tempCharge * chargeScale * elecScale;
              // Load the dipole in the global frame.
              out[t100] = dipole[0] * dipoleScale * elecScale;
              out[t010] = dipole[1] * dipoleScale * elecScale;
              out[t001] = dipole[2] * dipoleScale * elecScale;
              // Load the quadrupole in the global frame.
              out[t200] = quadrupole[0][0] * quadrupoleScale * elecScale;
              out[t020] = quadrupole[1][1] * quadrupoleScale * elecScale;
              out[t002] = quadrupole[2][2] * quadrupoleScale * elecScale;
              out[t110] = quadrupole[0][1] * quadrupoleScale * elecScale;
              out[t101] = quadrupole[0][2] * quadrupoleScale * elecScale;
              out[t011] = quadrupole[1][2] * quadrupoleScale * elecScale;

              tempCharge = esvMultipoleTautDot[0];
              // Load the dipole for rotation.
              tempDipole[0] = esvMultipoleTautDot[t100];
              tempDipole[1] = esvMultipoleTautDot[t010];
              tempDipole[2] = esvMultipoleTautDot[t001];
              // Load the quadrupole for rotation.
              tempQuadrupole[0][0] = esvMultipoleTautDot[t200];
              tempQuadrupole[1][1] = esvMultipoleTautDot[t020];
              tempQuadrupole[2][2] = esvMultipoleTautDot[t002];
              tempQuadrupole[0][1] = esvMultipoleTautDot[t110];
              tempQuadrupole[0][2] = esvMultipoleTautDot[t101];
              tempQuadrupole[1][2] = esvMultipoleTautDot[t011];
              tempQuadrupole[1][0] = esvMultipoleTautDot[t110];
              tempQuadrupole[2][0] = esvMultipoleTautDot[t101];
              tempQuadrupole[2][1] = esvMultipoleTautDot[t011];
              if (needsChiralInversion) {
                tempDipole[1] = -tempDipole[1];
                tempQuadrupole[0][1] = -tempQuadrupole[0][1];
                tempQuadrupole[1][0] = -tempQuadrupole[1][0];
                tempQuadrupole[1][2] = -tempQuadrupole[1][2];
                tempQuadrupole[2][1] = -tempQuadrupole[2][1];
              }
              rotateMultipole(rotmat, tempDipole, tempQuadrupole, dipole, quadrupole);
              out = tautomerMultipole[iSymm][ii];
              out[t000] = tempCharge * chargeScale * elecScale;
              // Load the dipole in the global frame.
              out[t100] = dipole[0] * dipoleScale * elecScale;
              out[t010] = dipole[1] * dipoleScale * elecScale;
              out[t001] = dipole[2] * dipoleScale * elecScale;
              // Load the quadrupole in the global frame.
              out[t200] = quadrupole[0][0] * quadrupoleScale * elecScale;
              out[t020] = quadrupole[1][1] * quadrupoleScale * elecScale;
              out[t002] = quadrupole[2][2] * quadrupoleScale * elecScale;
              out[t110] = quadrupole[0][1] * quadrupoleScale * elecScale;
              out[t101] = quadrupole[0][2] * quadrupoleScale * elecScale;
              out[t011] = quadrupole[1][2] * quadrupoleScale * elecScale;
            }
          } else {
            // No multipole rotation for isolating torque vs. non-torque pieces of the multipole
            // energy gradient.
            double[] out = globalMultipole[iSymm][ii];
            out[t000] = in[t000] * chargeScale * elecScale;
            out[t100] = in[t100] * dipoleScale * elecScale;
            out[t010] = in[t010] * dipoleScale * elecScale;
            out[t001] = in[t001] * dipoleScale * elecScale;
            out[t200] = in[t200] * quadrupoleScale * elecScale;
            out[t020] = in[t020] * quadrupoleScale * elecScale;
            out[t002] = in[t002] * quadrupoleScale * elecScale;
            out[t110] = in[t110] * quadrupoleScale * elecScale;
            out[t101] = in[t101] * quadrupoleScale * elecScale;
            out[t011] = in[t011] * quadrupoleScale * elecScale;
            /* For ESV atoms, also scale the Mdot multipole. */
            if (esvTerm && esvSystem.isTitrating(ii)) {
              double[] esvMultipoleTitrDot = new double[10];
              double titrationLambda = esvSystem.getTitrationLambda(ii);
              double tautomerLambda = esvSystem.getTautomerLambda(ii);
              esvMultipoleTitrDot = esvSystem.getTitrationUtils()
                  .getMultipoleTitrationDeriv(atom, titrationLambda,
                      tautomerLambda, esvMultipoleTitrDot);
              in = esvMultipoleTitrDot;
              out = titrationMultipole[iSymm][ii];
              out[t000] = in[t000] * chargeScale * elecScale;
              out[t100] = in[t100] * dipoleScale * elecScale;
              out[t010] = in[t010] * dipoleScale * elecScale;
              out[t001] = in[t001] * dipoleScale * elecScale;
              out[t200] = in[t200] * quadrupoleScale * elecScale;
              out[t020] = in[t020] * quadrupoleScale * elecScale;
              out[t002] = in[t002] * quadrupoleScale * elecScale;
              out[t110] = in[t110] * quadrupoleScale * elecScale;
              out[t101] = in[t101] * quadrupoleScale * elecScale;
              out[t011] = in[t011] * quadrupoleScale * elecScale;

              double[] esvMultipoleTautDot = new double[10];
              esvMultipoleTautDot = esvSystem.getTitrationUtils()
                  .getMultipoleTautomerDeriv(atom, titrationLambda,
                      tautomerLambda, esvMultipoleTautDot);
              in = esvMultipoleTautDot;
              out = tautomerMultipole[iSymm][ii];
              out[t000] = in[t000] * chargeScale * elecScale;
              out[t100] = in[t100] * dipoleScale * elecScale;
              out[t010] = in[t010] * dipoleScale * elecScale;
              out[t001] = in[t001] * dipoleScale * elecScale;
              out[t200] = in[t200] * quadrupoleScale * elecScale;
              out[t020] = in[t020] * quadrupoleScale * elecScale;
              out[t002] = in[t002] * quadrupoleScale * elecScale;
              out[t110] = in[t110] * quadrupoleScale * elecScale;
              out[t101] = in[t101] * quadrupoleScale * elecScale;
              out[t011] = in[t011] * quadrupoleScale * elecScale;
            }
          }

          // Load the polarizability.
          PolarizeType polarizeType = particleMeshEwald.getPolarizeType(ii);
          if (polarizeType != null) {
            polarizability[ii] = polarizeType.polarizability * polarizabilityScale * elecScale;
            if (esvTerm && esvSystem.isTitrating(ii) && (esvSystem.isTitratingHydrogen(ii)
                || esvSystem.isTitratingSulfur(ii))) {
              titrationPolarizability[ii] = 0.0;
              tautomerPolarizability[ii] = 0.0;
              double titrationLambda = esvSystem.getTitrationLambda(ii);
              double tautomerLambda = esvSystem.getTautomerLambda(ii);
              titrationPolarizability[ii] = esvSystem.getTitrationUtils()
                  .getPolarizabilityTitrationDeriv(atom, titrationLambda, tautomerLambda);
              tautomerPolarizability[ii] = esvSystem.getTitrationUtils()
                  .getPolarizabilityTautomerDeriv(atom, titrationLambda, tautomerLambda);
            }
            thole[ii] = polarizeType.thole;
            ipdamp[ii] = polarizeType.pdamp;
            if (!(ipdamp[ii] > 0.0)) {
              ipdamp[ii] = Double.POSITIVE_INFINITY;
            } else {
              ipdamp[ii] = 1.0 / ipdamp[ii];
            }
          } else {
            polarizability[ii] = 0.0;
            thole[ii] = 0.0;
            ipdamp[ii] = Double.POSITIVE_INFINITY;
          }
        }
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return IntegerSchedule.fixed();
    }
  }
}
