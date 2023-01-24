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
package ffx.potential.nonbonded.pme;

import static ffx.potential.parameters.MultipoleType.t000;
import static ffx.potential.parameters.MultipoleType.t001;
import static ffx.potential.parameters.MultipoleType.t002;
import static ffx.potential.parameters.MultipoleType.t003;
import static ffx.potential.parameters.MultipoleType.t010;
import static ffx.potential.parameters.MultipoleType.t011;
import static ffx.potential.parameters.MultipoleType.t012;
import static ffx.potential.parameters.MultipoleType.t020;
import static ffx.potential.parameters.MultipoleType.t021;
import static ffx.potential.parameters.MultipoleType.t030;
import static ffx.potential.parameters.MultipoleType.t100;
import static ffx.potential.parameters.MultipoleType.t101;
import static ffx.potential.parameters.MultipoleType.t102;
import static ffx.potential.parameters.MultipoleType.t110;
import static ffx.potential.parameters.MultipoleType.t111;
import static ffx.potential.parameters.MultipoleType.t120;
import static ffx.potential.parameters.MultipoleType.t200;
import static ffx.potential.parameters.MultipoleType.t201;
import static ffx.potential.parameters.MultipoleType.t210;
import static ffx.potential.parameters.MultipoleType.t300;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;
import ffx.crystal.Crystal;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.numerics.multipole.MultipoleTensor;
import ffx.potential.bonded.Atom;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.nonbonded.ReciprocalSpace;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Parallel evaluation of the PME reciprocal space energy and gradient.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ReciprocalEnergyRegion extends ParallelRegion {

  private static final Logger logger = Logger.getLogger(ReciprocalEnergyRegion.class.getName());

  private static final double SQRT_PI = sqrt(Math.PI);
  private static final int tensorCount = MultipoleTensor.tensorCount(3);
  private final double electric;
  private final double aewald1;
  private final double aewald2;
  private final double aewald3;
  private final double aewald4;
  private final double oneThird;
  private final double twoThirds;
  private final int maxThreads;
  private final SharedDouble inducedDipoleSelfEnergy;
  private final SharedDouble inducedDipoleRecipEnergy;
  private final PermanentReciprocalEnergyLoop[] permanentReciprocalEnergyLoop;
  private final InducedDipoleReciprocalEnergyLoop[] inducedDipoleReciprocalEnergyLoop;
  /** Dimensions of [nsymm][nAtoms][3] */
  public double[][][] inducedDipole;
  public double[][][] inducedDipoleCR;
  /** Convenience reference to inducedDipole[0] */
  private double[][] ind;
  /** Convenience reference to inducedDipoleCR[0] */
  private double[][] indCR;
  private double nfftX, nfftY, nfftZ;
  private double[][] multipole;
  private double[][] fracMultipole;
  private double[][] fracMultipolePhi;
  private double[][] fracInducedDipolePhi;
  private double[][] fracInducedDipolePhiCR;
  private double permanentSelfEnergy;
  private double permanentReciprocalEnergy;
  /** An ordered array of atoms in the system. */
  private Atom[] atoms;
  /** Unit cell and spacegroup information. */
  private Crystal crystal;
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
  /** Dimensions of [nsymm][nAtoms][10] */
  private double[][][] globalMultipole;
  private double[][][] globalFracMultipole;
  private double[][][] titrationMultipole;
  private double[][][] tautomerMultipole;
  private ExtendedSystem extendedSystem = null;
  private boolean esvTerm = false;
  private double[][] cartMultipolePhi;
  private double[][] cartesianDipolePhi;
  private double[][] cartesianDipolePhiCR;
  /** Reciprocal space instance. */
  private ReciprocalSpace reciprocalSpace;
  private Polarization polarization;
  /** Atomic Gradient array. */
  private AtomicDoubleArray3D grad;
  /** Atomic Torque array. */
  private AtomicDoubleArray3D torque;
  /** Partial derivative of the gradient with respect to Lambda. */
  private AtomicDoubleArray3D lambdaGrad;
  /** Partial derivative of the torque with respect to Lambda. */
  private AtomicDoubleArray3D lambdaTorque;
  /** If true, compute coordinate gradient. */
  private boolean gradient;
  /**
   * If lambdaTerm is true, some ligand atom interactions with the environment are being turned
   * on/off.
   */
  private boolean lambdaTerm;
  /** Partial derivative with respect to Lambda. */
  private SharedDouble shareddEdLambda;
  /** Second partial derivative with respect to Lambda. */
  private SharedDouble sharedd2EdLambda2;
  // Alchemical control
  private double permanentScale;
  private double dlPowPerm;
  private double d2lPowPerm;
  private double polarizationScale;
  private double dlPowPol;
  private double d2lPowPol;
  private double dEdLSign;

  public ReciprocalEnergyRegion(int nt, double aewald, double electric) {
    permanentReciprocalEnergyLoop = new PermanentReciprocalEnergyLoop[nt];
    inducedDipoleReciprocalEnergyLoop = new InducedDipoleReciprocalEnergyLoop[nt];
    inducedDipoleSelfEnergy = new SharedDouble();
    inducedDipoleRecipEnergy = new SharedDouble();
    maxThreads = nt;

    this.electric = electric;
    aewald1 = -electric * aewald / SQRT_PI;
    aewald2 = 2.0 * aewald * aewald;
    aewald3 = -2.0 / 3.0 * electric * aewald * aewald * aewald / SQRT_PI;
    aewald4 = -2.0 * aewald3;
    oneThird = 1.0 / 3.0;
    twoThirds = 2.0 / 3.0;
  }

  /**
   * Execute the ReciprocalEnergyRegion with the passed ParallelTeam.
   *
   * @param parallelTeam The ParallelTeam instance to execute with.
   */
  public void executeWith(ParallelTeam parallelTeam) {
    try {
      parallelTeam.execute(this);
    } catch (Exception e) {
      String message = " Exception computing the electrostatic energy.\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  @Override
  public void finish() {
    /*
     The permanent multipole self energy contributions are large
     enough that rounding differences that result from threads
     finishing in different orders removes deterministic behavior.
    */
    permanentSelfEnergy = 0.0;
    permanentReciprocalEnergy = 0.0;
    for (int i = 0; i < maxThreads; i++) {
      permanentSelfEnergy += permanentReciprocalEnergyLoop[i].eSelf;
      permanentReciprocalEnergy += permanentReciprocalEnergyLoop[i].eRecip;
    }
  }

  public double getInducedDipoleReciprocalEnergy() {
    return inducedDipoleRecipEnergy.get();
  }

  public double getInducedDipoleSelfEnergy() {
    return inducedDipoleSelfEnergy.get();
  }

  public double getPermanentReciprocalEnergy() {
    return permanentReciprocalEnergy;
  }

  public double getPermanentSelfEnergy() {
    return permanentSelfEnergy;
  }

  public void init(
      Atom[] atoms,
      Crystal crystal,
      boolean gradient,
      boolean lambdaTerm,
      boolean esvTerm,
      boolean[] use,
      double[][][] globalMultipole,
      double[][][] globalFracMultipole,
      double[][][] titrationMultipole,
      double[][][] tautomerMultipole,
      double[][] cartMultipolePhi,
      double[][] fracMultipolePhi,
      Polarization polarization,
      double[][][] inducedDipole,
      double[][][] inducedDipoleCR,
      double[][] cartesianDipolePhi,
      double[][] cartesianDipolePhiCR,
      double[][] fracInducedDipolePhi,
      double[][] fracInducedDipolePhiCR,
      ReciprocalSpace reciprocalSpace,
      AlchemicalParameters alchemicalParameters,
      ExtendedSystem extendedSystem,
      AtomicDoubleArray3D grad,
      AtomicDoubleArray3D torque,
      AtomicDoubleArray3D lambdaGrad,
      AtomicDoubleArray3D lambdaTorque,
      SharedDouble shareddEdLambda,
      SharedDouble sharedd2EdLambda2
      ) {
    this.atoms = atoms;
    this.crystal = crystal;
    this.gradient = gradient;
    this.lambdaTerm = lambdaTerm;
    this.esvTerm = esvTerm;
    this.use = use;
    // Permanent
    this.globalMultipole = globalMultipole;
    this.globalFracMultipole = globalFracMultipole;
    this.titrationMultipole = titrationMultipole;
    this.tautomerMultipole = tautomerMultipole;
    this.cartMultipolePhi = cartMultipolePhi;
    this.fracMultipolePhi = fracMultipolePhi;
    // Polarization
    this.polarization = polarization;
    this.inducedDipole = inducedDipole;
    this.inducedDipoleCR = inducedDipoleCR;
    this.cartesianDipolePhi = cartesianDipolePhi;
    this.cartesianDipolePhiCR = cartesianDipolePhiCR;
    this.fracInducedDipolePhi = fracInducedDipolePhi;
    this.fracInducedDipolePhiCR = fracInducedDipolePhiCR;
    this.reciprocalSpace = reciprocalSpace;
    this.extendedSystem = extendedSystem;
    // Alchemical control
    this.permanentScale = alchemicalParameters.permanentScale;
    this.dlPowPerm = alchemicalParameters.dlPowPerm;
    this.d2lPowPerm = alchemicalParameters.d2lPowPerm;
    this.polarizationScale = alchemicalParameters.polarizationScale;
    this.dlPowPol = alchemicalParameters.dlPowPol;
    this.d2lPowPol = alchemicalParameters.d2lPowPol;
    this.dEdLSign = alchemicalParameters.dEdLSign;
    // Output
    this.grad = grad;
    this.torque = torque;
    this.lambdaGrad = lambdaGrad;
    this.lambdaTorque = lambdaTorque;
    this.shareddEdLambda = shareddEdLambda;
    this.sharedd2EdLambda2 = sharedd2EdLambda2;
  }

  @Override
  public void run() throws Exception {
    int threadIndex = getThreadIndex();
    if (permanentReciprocalEnergyLoop[threadIndex] == null) {
      permanentReciprocalEnergyLoop[threadIndex] = new PermanentReciprocalEnergyLoop();
      inducedDipoleReciprocalEnergyLoop[threadIndex] = new InducedDipoleReciprocalEnergyLoop();
    }
    try {
      int nAtoms = atoms.length;
      execute(0, nAtoms - 1, permanentReciprocalEnergyLoop[threadIndex]);
      if (polarization != Polarization.NONE) {
        execute(0, nAtoms - 1, inducedDipoleReciprocalEnergyLoop[threadIndex]);
      }
    } catch (Exception e) {
      String message =
          "Fatal exception computing the real space field in thread " + threadIndex + "\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  @Override
  public void start() {
    multipole = globalMultipole[0];
    fracMultipole = globalFracMultipole[0];
    ind = inducedDipole[0];
    indCR = inducedDipoleCR[0];
    inducedDipoleSelfEnergy.set(0.0);
    inducedDipoleRecipEnergy.set(0.0);
    nfftX = reciprocalSpace.getXDim();
    nfftY = reciprocalSpace.getYDim();
    nfftZ = reciprocalSpace.getZDim();
  }

  private class PermanentReciprocalEnergyLoop extends IntegerForLoop {

    int threadID;
    double eSelf;
    double eRecip;

    @Override
    public void finish() {
      eSelf *= permanentScale;
      eRecip *= permanentScale * 0.5 * electric;
    }

    @Override
    public void run(int lb, int ub) throws Exception {

      // Permanent multipole self energy and gradient.
      for (int i = lb; i <= ub; i++) {
        if (use[i]) {
          double[] in = globalMultipole[0][i];
          double cii = in[t000] * in[t000];
          double dii = in[t100] * in[t100] + in[t010] * in[t010] + in[t001] * in[t001];
          double qii = in[t200] * in[t200] + in[t020] * in[t020] + in[t002] * in[t002]
              + 2.0 * (in[t110] * in[t110] + in[t101] * in[t101] + in[t011] * in[t011]);
          eSelf += aewald1 * (cii + aewald2 * (dii / 3.0 + 2.0 * aewald2 * qii / 45.0));

          if (esvTerm && extendedSystem.isTitrating(i)) {
            double[] indot = titrationMultipole[0][i];
            double ciidot = indot[t000] * in[t000];
            double diidot = indot[t100] * in[t100] + indot[t010] * in[t010] + indot[t001] * in[t001];
            double qiidot = indot[t200] * in[t200] + indot[t020] * in[t020] + indot[t002] * in[t002]
                + 2.0 * (indot[t110] * in[t110] + indot[t101] * in[t101] + indot[t011] * in[t011]);
            double eSelfTitrDot =
                2.0 * aewald1 * (ciidot + aewald2 * (diidot / 3.0 + 2.0 * aewald2 * qiidot / 45.0));
            double eSelfTautDot = 0.0;
            if (extendedSystem.isTautomerizing(i)) {
              indot = tautomerMultipole[0][i];
              ciidot = indot[t000] * in[t000];
              diidot = indot[t100] * in[t100] + indot[t010] * in[t010] + indot[t001] * in[t001];
              qiidot = indot[t200] * in[t200] + indot[t020] * in[t020] + indot[t002] * in[t002]
                  + 2.0 * (indot[t110] * in[t110] + indot[t101] * in[t101] + indot[t011] * in[t011]);
              eSelfTautDot = 2.0 * aewald1 * (ciidot + aewald2 * (diidot / 3.0
                  + 2.0 * aewald2 * qiidot / 45.0));
            }
            double factor = permanentScale;
            extendedSystem.addPermElecDeriv(i, eSelfTitrDot * factor, eSelfTautDot * factor);
          }
        }
      }
      if (lambdaTerm) {
        shareddEdLambda.addAndGet(eSelf * dlPowPerm * dEdLSign);
        sharedd2EdLambda2.addAndGet(eSelf * d2lPowPerm * dEdLSign);
      }

      // Permanent multipole reciprocal space energy and gradient.
      final double[][] recip = crystal.getUnitCell().A;

      double dUdL = 0.0;
      double d2UdL2 = 0.0;
      for (int i = lb; i <= ub; i++) {
        if (use[i]) {
          final double[] phi = cartMultipolePhi[i];
          final double[] mpole = multipole[i];
          final double[] fmpole = fracMultipole[i];
          double e = mpole[t000] * phi[t000]
              + mpole[t100] * phi[t100] + mpole[t010] * phi[t010] + mpole[t001] * phi[t001]
              + oneThird * (mpole[t200] * phi[t200] + mpole[t020] * phi[t020]
              + mpole[t002] * phi[t002]
              + 2.0 * (mpole[t110] * phi[t110] + mpole[t101] * phi[t101] + mpole[t011] * phi[t011]));
          eRecip += e;
          if (esvTerm && extendedSystem.isTitrating(i)) {
            double[] mpoleDot = titrationMultipole[0][i];
            double edotTitr = mpoleDot[t000] * phi[t000]
                + mpoleDot[t100] * phi[t100] + mpoleDot[t010] * phi[t010]
                + mpoleDot[t001] * phi[t001]
                + oneThird * (mpoleDot[t200] * phi[t200] + mpoleDot[t020] * phi[t020]
                + mpoleDot[t002] * phi[t002]
                + 2.0 * (mpoleDot[t110] * phi[t110] + mpoleDot[t101] * phi[t101]
                + mpoleDot[t011] * phi[t011]));
            double edotTaut = 0.0;
            if (extendedSystem.isTautomerizing(i)) {
              mpoleDot = tautomerMultipole[0][i];
              edotTaut = mpoleDot[t000] * phi[t000]
                  + mpoleDot[t100] * phi[t100] + mpoleDot[t010] * phi[t010]
                  + mpoleDot[t001] * phi[t001]
                  + oneThird * (mpoleDot[t200] * phi[t200] + mpoleDot[t020] * phi[t020]
                  + mpoleDot[t002] * phi[t002]
                  + 2.0 * (mpoleDot[t110] * phi[t110] + mpoleDot[t101] * phi[t101]
                  + mpoleDot[t011] * phi[t011]));
            }
            double factor = permanentScale * electric;
            extendedSystem.addPermElecDeriv(i, edotTitr * factor, edotTaut * factor);
          }
          if (gradient || lambdaTerm) {
            final double[] fPhi = fracMultipolePhi[i];
            double gx = fmpole[t000] * fPhi[t100]
                + fmpole[t100] * fPhi[t200] + fmpole[t010] * fPhi[t110] + fmpole[t001] * fPhi[t101]
                + fmpole[t200] * fPhi[t300] + fmpole[t020] * fPhi[t120] + fmpole[t002] * fPhi[t102]
                + fmpole[t110] * fPhi[t210] + fmpole[t101] * fPhi[t201] + fmpole[t011] * fPhi[t111];
            double gy = fmpole[t000] * fPhi[t010]
                + fmpole[t100] * fPhi[t110] + fmpole[t010] * fPhi[t020] + fmpole[t001] * fPhi[t011]
                + fmpole[t200] * fPhi[t210] + fmpole[t020] * fPhi[t030] + fmpole[t002] * fPhi[t012]
                + fmpole[t110] * fPhi[t120] + fmpole[t101] * fPhi[t111] + fmpole[t011] * fPhi[t021];
            double gz = fmpole[t000] * fPhi[t001]
                + fmpole[t100] * fPhi[t101] + fmpole[t010] * fPhi[t011] + fmpole[t001] * fPhi[t002]
                + fmpole[t200] * fPhi[t201] + fmpole[t020] * fPhi[t021] + fmpole[t002] * fPhi[t003]
                + fmpole[t110] * fPhi[t111] + fmpole[t101] * fPhi[t102] + fmpole[t011] * fPhi[t012];
            gx *= nfftX;
            gy *= nfftY;
            gz *= nfftZ;
            final double dfx = recip[0][0] * gx + recip[0][1] * gy + recip[0][2] * gz;
            final double dfy = recip[1][0] * gx + recip[1][1] * gy + recip[1][2] * gz;
            final double dfz = recip[2][0] * gx + recip[2][1] * gy + recip[2][2] * gz;
            // Compute dipole torques
            double tqx = -mpole[t010] * phi[t001] + mpole[t001] * phi[t010];
            double tqy = -mpole[t001] * phi[t100] + mpole[t100] * phi[t001];
            double tqz = -mpole[t100] * phi[t010] + mpole[t010] * phi[t100];
            // Compute quadrupole torques
            tqx -= twoThirds * (
                mpole[t110] * phi[t101] + mpole[t020] * phi[t011] + mpole[t011] * phi[t002]
                    - mpole[t101] * phi[t110] - mpole[t011] * phi[t020] - mpole[t002] * phi[t011]);
            tqy -= twoThirds * (
                mpole[t101] * phi[t200] + mpole[t011] * phi[t110] + mpole[t002] * phi[t101]
                    - mpole[t200] * phi[t101] - mpole[t110] * phi[t011] - mpole[t101] * phi[t002]);
            tqz -= twoThirds * (
                mpole[t200] * phi[t110] + mpole[t110] * phi[t020] + mpole[t101] * phi[t011]
                    - mpole[t110] * phi[t200] - mpole[t020] * phi[t110] - mpole[t011] * phi[t101]);
            if (gradient) {
              double factor = permanentScale * electric;
              grad.add(threadID, i, factor * dfx, factor * dfy, factor * dfz);
              torque.add(threadID, i, factor * tqx, factor * tqy, factor * tqz);
            }
            if (lambdaTerm) {
              dUdL += dEdLSign * dlPowPerm * e;
              d2UdL2 += dEdLSign * d2lPowPerm * e;
              double factor = dEdLSign * dlPowPerm * electric;
              lambdaGrad.add(threadID, i, factor * dfx, factor * dfy, factor * dfz);
              lambdaTorque.add(threadID, i, factor * tqx, factor * tqy, factor * tqz);
            }
          }
        }
      }

      if (lambdaTerm) {
        shareddEdLambda.addAndGet(0.5 * dUdL * electric);
        sharedd2EdLambda2.addAndGet(0.5 * d2UdL2 * electric);
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return IntegerSchedule.fixed();
    }

    @Override
    public void start() {
      eSelf = 0.0;
      eRecip = 0.0;
      threadID = getThreadIndex();
    }
  }

  private class InducedDipoleReciprocalEnergyLoop extends IntegerForLoop {

    private final double[] sfPhi = new double[tensorCount];
    private final double[] sPhi = new double[tensorCount];
    private double eSelf;
    private double eRecip;
    private int threadID;

    @Override
    public void finish() {
      inducedDipoleSelfEnergy.addAndGet(polarizationScale * eSelf);
      inducedDipoleRecipEnergy.addAndGet(polarizationScale * eRecip);
    }

    @Override
    public void run(int lb, int ub) throws Exception {
      // Induced dipole self energy and gradient.
      for (int i = lb; i <= ub; i++) {
        if (use[i]) {
          final double[] indi = ind[i];
          final double[] multipolei = multipole[i];
          final double dix = multipolei[t100];
          final double diy = multipolei[t010];
          final double diz = multipolei[t001];
          final double dii = indi[0] * dix + indi[1] * diy + indi[2] * diz;
          eSelf += aewald3 * dii;
        }
      }
      if (lambdaTerm) {
        shareddEdLambda.addAndGet(dEdLSign * dlPowPol * eSelf);
        sharedd2EdLambda2.addAndGet(dEdLSign * d2lPowPol * eSelf);
      }
      if (gradient || esvTerm) {
        for (int i = lb; i <= ub; i++) {
          if (use[i]) {
            final double[] indi = ind[i];
            final double[] indpi = indCR[i];
            final double[] multipolei = multipole[i];
            final double dix = multipolei[t100];
            final double diy = multipolei[t010];
            final double diz = multipolei[t001];
            final double uix = 0.5 * (indi[0] + indpi[0]);
            final double uiy = 0.5 * (indi[1] + indpi[1]);
            final double uiz = 0.5 * (indi[2] + indpi[2]);
            if (gradient) {
              final double tix = aewald4 * (diy * uiz - diz * uiy);
              final double tiy = aewald4 * (diz * uix - dix * uiz);
              final double tiz = aewald4 * (dix * uiy - diy * uix);
              torque.add(threadID, i, polarizationScale * tix, polarizationScale * tiy,
                  polarizationScale * tiz);
              if (lambdaTerm) {
                double factor = dEdLSign * dlPowPol;
                lambdaTorque.add(threadID, i, factor * tix, factor * tiy, factor * tiz);
              }
            }
            if (esvTerm && extendedSystem.isTitrating(i)) {
              double[] titrDot = titrationMultipole[0][i];
              double indiDotMdot = uix * titrDot[t100] + uiy * titrDot[t010] + uiz * titrDot[t001];
              double dIndSelfTitrdLi = polarizationScale * 2.0 * aewald3 * indiDotMdot;
              double dIndSelfTautdLi = 0.0;
              if (extendedSystem.isTautomerizing(i)) {
                double[] tautDot = tautomerMultipole[0][i];
                indiDotMdot = uix * tautDot[t100] + uiy * tautDot[t010] + uiz * tautDot[t001];
                dIndSelfTautdLi = polarizationScale * 2.0 * aewald3 * indiDotMdot;
              }
              extendedSystem.addIndElecDeriv(i, dIndSelfTitrdLi, dIndSelfTautdLi);
            }
          }
        }
      }

      // Induced dipole reciprocal space energy and gradient.
      for (int i = lb; i <= ub; i++) {
        double[] fracInd = new double[3];
        double[] fracIndCR = new double[3];
        if (use[i]) {
          final double[] fPhi = fracMultipolePhi[i];
          reciprocalSpace.cartToFracInducedDipole(inducedDipole[0][i], inducedDipoleCR[0][i], fracInd, fracIndCR);
          final double indx = fracInd[0];
          final double indy = fracInd[1];
          final double indz = fracInd[2];
          eRecip += indx * fPhi[t100] + indy * fPhi[t010] + indz * fPhi[t001];
          if (gradient || esvTerm) {
            final double[] iPhi = cartesianDipolePhi[i];
            final double[] iCRPhi = cartesianDipolePhiCR[i];
            final double[] fiPhi = fracInducedDipolePhi[i];
            final double[] fiCRPhi = fracInducedDipolePhiCR[i];
            final double[] mpolei = multipole[i];
            final double[] fmpolei = fracMultipole[i];
            final double inpx = fracIndCR[0];
            final double inpy = fracIndCR[1];
            final double inpz = fracIndCR[2];
            final double insx = indx + inpx;
            final double insy = indy + inpy;
            final double insz = indz + inpz;
            for (int t = 0; t < tensorCount; t++) {
              sPhi[t] = 0.5 * (iPhi[t] + iCRPhi[t]);
              sfPhi[t] = fiPhi[t] + fiCRPhi[t];
            }
            double gx = insx * fPhi[t200] + insy * fPhi[t110] + insz * fPhi[t101];
            double gy = insx * fPhi[t110] + insy * fPhi[t020] + insz * fPhi[t011];
            double gz = insx * fPhi[t101] + insy * fPhi[t011] + insz * fPhi[t002];
            if (polarization == Polarization.MUTUAL) {
              gx += indx * fiCRPhi[t200] + inpx * fiPhi[t200] + indy * fiCRPhi[t110]
                  + inpy * fiPhi[t110] + indz * fiCRPhi[t101] + inpz * fiPhi[t101];
              gy += indx * fiCRPhi[t110] + inpx * fiPhi[t110] + indy * fiCRPhi[t020]
                  + inpy * fiPhi[t020] + indz * fiCRPhi[t011] + inpz * fiPhi[t011];
              gz += indx * fiCRPhi[t101] + inpx * fiPhi[t101] + indy * fiCRPhi[t011]
                  + inpy * fiPhi[t011] + indz * fiCRPhi[t002] + inpz * fiPhi[t002];
            }
            gx += fmpolei[t000] * sfPhi[t100] + fmpolei[t100] * sfPhi[t200]
                + fmpolei[t010] * sfPhi[t110] + fmpolei[t001] * sfPhi[t101]
                + fmpolei[t200] * sfPhi[t300] + fmpolei[t020] * sfPhi[t120]
                + fmpolei[t002] * sfPhi[t102] + fmpolei[t110] * sfPhi[t210]
                + fmpolei[t101] * sfPhi[t201] + fmpolei[t011] * sfPhi[t111];
            gy += fmpolei[t000] * sfPhi[t010] + fmpolei[t100] * sfPhi[t110]
                + fmpolei[t010] * sfPhi[t020] + fmpolei[t001] * sfPhi[t011]
                + fmpolei[t200] * sfPhi[t210] + fmpolei[t020] * sfPhi[t030]
                + fmpolei[t002] * sfPhi[t012] + fmpolei[t110] * sfPhi[t120]
                + fmpolei[t101] * sfPhi[t111] + fmpolei[t011] * sfPhi[t021];
            gz += fmpolei[t000] * sfPhi[t001] + fmpolei[t100] * sfPhi[t101]
                + fmpolei[t010] * sfPhi[t011] + fmpolei[t001] * sfPhi[t002]
                + fmpolei[t200] * sfPhi[t201] + fmpolei[t020] * sfPhi[t021]
                + fmpolei[t002] * sfPhi[t003] + fmpolei[t110] * sfPhi[t111]
                + fmpolei[t101] * sfPhi[t102] + fmpolei[t011] * sfPhi[t012];
            gx *= nfftX;
            gy *= nfftY;
            gz *= nfftZ;
            double[][] recip = crystal.getUnitCell().A;
            double dfx = recip[0][0] * gx + recip[0][1] * gy + recip[0][2] * gz;
            double dfy = recip[1][0] * gx + recip[1][1] * gy + recip[1][2] * gz;
            double dfz = recip[2][0] * gx + recip[2][1] * gy + recip[2][2] * gz;
            dfx *= 0.5 * electric;
            dfy *= 0.5 * electric;
            dfz *= 0.5 * electric;
            // Compute dipole torques
            double tqx = -mpolei[t010] * sPhi[t001] + mpolei[t001] * sPhi[t010];
            double tqy = -mpolei[t001] * sPhi[t100] + mpolei[t100] * sPhi[t001];
            double tqz = -mpolei[t100] * sPhi[t010] + mpolei[t010] * sPhi[t100];
            // Compute quadrupole torques
            tqx -= twoThirds * (
                mpolei[t110] * sPhi[t101] + mpolei[t020] * sPhi[t011] + mpolei[t011] * sPhi[t002]
                    - mpolei[t101] * sPhi[t110] - mpolei[t011] * sPhi[t020]
                    - mpolei[t002] * sPhi[t011]);
            tqy -= twoThirds * (
                mpolei[t101] * sPhi[t200] + mpolei[t011] * sPhi[t110] + mpolei[t002] * sPhi[t101]
                    - mpolei[t200] * sPhi[t101] - mpolei[t110] * sPhi[t011]
                    - mpolei[t101] * sPhi[t002]);
            tqz -= twoThirds * (
                mpolei[t200] * sPhi[t110] + mpolei[t110] * sPhi[t020] + mpolei[t101] * sPhi[t011]
                    - mpolei[t110] * sPhi[t200] - mpolei[t020] * sPhi[t110]
                    - mpolei[t011] * sPhi[t101]);
            tqx *= electric;
            tqy *= electric;
            tqz *= electric;
            grad.add(threadID, i, polarizationScale * dfx, polarizationScale * dfy,
                polarizationScale * dfz);
            torque.add(threadID, i, polarizationScale * tqx, polarizationScale * tqy,
                polarizationScale * tqz);
            if (lambdaTerm) {
              double factor = dEdLSign * dlPowPol;
              lambdaGrad.add(threadID, i, factor * dfx, factor * dfy, factor * dfz);
              lambdaTorque.add(threadID, i, factor * tqx, factor * tqy, factor * tqz);
            }
            if (esvTerm && extendedSystem.isTitrating(i)) {
              final double[] mTitrDot = titrationMultipole[0][i];
              double eTitrDot = mTitrDot[t000] * sPhi[t000]
                  + mTitrDot[t100] * sPhi[t100] + mTitrDot[t010] * sPhi[t010]
                  + mTitrDot[t001] * sPhi[t001]
                  + oneThird * (mTitrDot[t200] * sPhi[t200] + mTitrDot[t020] * sPhi[t020]
                  + mTitrDot[t002] * sPhi[t002]
                  + 2.0 * (mTitrDot[t110] * sPhi[t110] + mTitrDot[t101] * sPhi[t101]
                  + mTitrDot[t011] * sPhi[t011]));
              double eTautDot = 0.0;
              if (extendedSystem.isTautomerizing(i)) {
                final double[] mTautDot = tautomerMultipole[0][i];
                eTautDot = mTautDot[t000] * sPhi[t000]
                    + mTautDot[t100] * sPhi[t100] + mTautDot[t010] * sPhi[t010]
                    + mTautDot[t001] * sPhi[t001]
                    + oneThird * (mTautDot[t200] * sPhi[t200] + mTautDot[t020] * sPhi[t020]
                    + mTautDot[t002] * sPhi[t002]
                    + 2.0 * (mTautDot[t110] * sPhi[t110]
                    + mTautDot[t101] * sPhi[t101] + mTautDot[t011] * sPhi[t011]));
              }
              double factor = polarizationScale * electric;
              extendedSystem.addIndElecDeriv(i, eTitrDot * factor, eTautDot * factor);
            }
          }
        }
      }
      eRecip *= 0.5 * electric;
      if (lambdaTerm) {
        shareddEdLambda.addAndGet(dEdLSign * dlPowPol * eRecip);
        sharedd2EdLambda2.addAndGet(dEdLSign * d2lPowPol * eRecip);
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return IntegerSchedule.fixed();
    }

    @Override
    public void start() {
      eSelf = 0.0;
      eRecip = 0.0;
      threadID = getThreadIndex();
    }
  }
}
