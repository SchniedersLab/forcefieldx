// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import ffx.potential.utils.EnergyException;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;

/**
 * Parallel conversion of torques into forces, and then reduce them.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ReduceRegion extends ParallelRegion {

  private static final Logger logger = Logger.getLogger(ReduceRegion.class.getName());

  /**
   * If set to false, multipoles are fixed in their local frame and torques are zero, which is
   * useful for narrowing down discrepancies between analytic and finite-difference
   * derivatives(default is true).
   */
  private final boolean rotateMultipoles;

  private final TorqueLoop[] torqueLoop;
  private final ReduceLoop[] reduceLoop;
  /**
   * If lambdaTerm is true, some ligand atom interactions with the environment are being turned
   * on/off.
   */
  private boolean lambdaTerm;
  /**
   * If true, compute coordinate gradient.
   */
  private boolean gradient;
  /**
   * An ordered array of atoms in the system.
   */
  private Atom[] atoms;
  /**
   * Dimensions of [nsymm][xyz][nAtoms].
   */
  private double[][][] coordinates;
  /**
   * Multipole frame definition.
   */
  private MultipoleFrameDefinition[] frame;
  /**
   * Multipole frame defining atoms.
   */
  private int[][] axisAtom;
  /**
   * Atomic Gradient array.
   */
  private AtomicDoubleArray3D grad;
  /**
   * Atomic Torque array.
   */
  private AtomicDoubleArray3D torque;
  /**
   * Partial derivative of the gradient with respect to Lambda.
   */
  private AtomicDoubleArray3D lambdaGrad;
  /**
   * Partial derivative of the torque with respect to Lambda.
   */
  private AtomicDoubleArray3D lambdaTorque;

  public ReduceRegion(int threadCount, ForceField forceField) {
    torqueLoop = new TorqueLoop[threadCount];
    reduceLoop = new ReduceLoop[threadCount];
    rotateMultipoles = forceField.getBoolean("ROTATE_MULTIPOLES", true);
  }

  /**
   * Execute the ReduceRegion with the passed ParallelTeam.
   *
   * @param parallelTeam The ParallelTeam instance to execute with.
   */
  public void executeWith(ParallelTeam parallelTeam) {
    try {
      parallelTeam.execute(this);
    } catch (Exception e) {
      String message = "Exception calculating torques.";
      logger.log(Level.SEVERE, message, e);
    }
  }

  public void init(
      boolean lambdaTerm,
      boolean gradient,
      Atom[] atoms,
      double[][][] coordinates,
      MultipoleFrameDefinition[] frame,
      int[][] axisAtom,
      AtomicDoubleArray3D grad,
      AtomicDoubleArray3D torque,
      AtomicDoubleArray3D lambdaGrad,
      AtomicDoubleArray3D lambdaTorque) {
    this.lambdaTerm = lambdaTerm;
    this.gradient = gradient;
    this.atoms = atoms;
    this.coordinates = coordinates;
    this.frame = frame;
    this.axisAtom = axisAtom;
    this.grad = grad;
    this.torque = torque;
    this.lambdaGrad = lambdaGrad;
    this.lambdaTorque = lambdaTorque;
  }

  @Override
  public void run() throws EnergyException {
    int nAtoms = atoms.length;
    try {
      int threadIndex = getThreadIndex();
      if (torqueLoop[threadIndex] == null) {
        torqueLoop[threadIndex] = new TorqueLoop();
        reduceLoop[threadIndex] = new ReduceLoop();
      }
      if (rotateMultipoles) {
        execute(0, nAtoms - 1, torqueLoop[threadIndex]);
      }
      execute(0, nAtoms - 1, reduceLoop[threadIndex]);
    } catch (Exception e) {
      if (e instanceof EnergyException) {
        throw (EnergyException) e;
      }
      String message = "Fatal exception computing torque in thread " + getThreadIndex() + "\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  private class TorqueLoop extends IntegerForLoop {

    private Torque torques;
    private int threadID;
    double[] trq = new double[3];
    double[][] g = new double[4][3];
    int[] frameIndex = new int[4];

    TorqueLoop() {
      torques = new Torque();
    }

    @Override
    public void run(int lb, int ub) throws EnergyException {
      if (gradient) {
        torque.reduce(lb, ub);
        for (int i = lb; i <= ub; i++) {
          // Gradients from torques will exist if the frameIndex is not -1.
          Arrays.fill(frameIndex, -1);
          trq[0] = torque.getX(i);
          trq[1] = torque.getY(i);
          trq[2] = torque.getZ(i);

          // Check for undefined torques.
          if (isNaN(trq[0]) || isInfinite(trq[0])
              || isNaN(trq[1]) || isInfinite(trq[1])
              || isNaN(trq[2]) || isInfinite(trq[2])) {
            Atom a = atoms[i];
            throw new EnergyException(
                format(" Undefined torque (%8.3f,%8.3f,%8.3f) for atom %s.", trq[0], trq[1], trq[2], a));
          }

          torques.torque(i, 0, trq, frameIndex, g);
          for (int j = 0; j < 4; j++) {
            int index = frameIndex[j];
            if (index >= 0) {
              double[] gj = g[j];

              // Check for undefined torques.
              if (isNaN(gj[0]) || isInfinite(gj[0])
                  || isNaN(gj[1]) || isInfinite(gj[1])
                  || isNaN(gj[2]) || isInfinite(gj[2])) {
                Atom ai = atoms[i];
                Atom aj = atoms[index];
                throw new EnergyException(
                    format(" Undefined gradient (%8.3f,%8.3f,%8.3f)\n For atom: %s\n From torque of atom %s",
                        gj[0], gj[1], gj[2], aj, ai));
              }


              grad.add(threadID, index, gj[0], gj[1], gj[2]);
            }
          }
        }
      }
      if (lambdaTerm) {
        lambdaTorque.reduce(lb, ub);
        for (int i = lb; i <= ub; i++) {
          // Gradients from torques will exist if the frameIndex is not -1.
          Arrays.fill(frameIndex, -1);
          trq[0] = lambdaTorque.getX(i);
          trq[1] = lambdaTorque.getY(i);
          trq[2] = lambdaTorque.getZ(i);
          torques.torque(i, 0, trq, frameIndex, g);
          for (int j = 0; j < 4; j++) {
            int index = frameIndex[j];
            if (index >= 0) {
              double[] gj = g[j];
              lambdaGrad.add(threadID, index, gj[0], gj[1], gj[2]);
            }
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
      threadID = getThreadIndex();
      torques.init(axisAtom, frame, coordinates);
    }
  }

  private class ReduceLoop extends IntegerForLoop {

    @Override
    public void run(int lb, int ub) throws EnergyException {
      if (gradient) {
        grad.reduce(lb, ub);
        for (int i = lb; i <= ub; i++) {
          Atom ai = atoms[i];
          ai.addToXYZGradient(grad.getX(i), grad.getY(i), grad.getZ(i));
        }
      }
      if (lambdaTerm) {
        lambdaGrad.reduce(lb, ub);
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return IntegerSchedule.fixed();
    }
  }
}
