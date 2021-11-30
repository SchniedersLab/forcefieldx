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
package ffx.potential.nonbonded.implicit;

import static java.lang.String.format;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import ffx.numerics.atomic.AtomicDoubleArray;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Parallel initialization of accumulation arrays for Generalized Kirkwood.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class InitializationRegion extends ParallelRegion {

  private static final Logger logger = Logger.getLogger(InitializationRegion.class.getName());

  /**
   * Reference to the GK instance so that parameters can be updated.
   */
  private GeneralizedKirkwood generalizedKirkwood;
  /**
   * Initialization loops.
   */
  private final InitializationLoop[] initializationLoop;
  /** Array of atoms. */
  private Atom[] atoms;
  /** True if GK is being turned on/off. */
  private boolean lambdaTerm;
  /** Atomic Gradient array. */
  private AtomicDoubleArray3D grad;
  /** Atomic Torque array. */
  private AtomicDoubleArray3D torque;
  /** Shared array for computation of Born radii gradient. */
  private AtomicDoubleArray sharedBornGrad;

  public InitializationRegion(int maxThreads) {
    initializationLoop = new InitializationLoop[maxThreads];
  }

  /**
   * Execute the InitializationRegion with the passed ParallelTeam.
   *
   * @param parallelTeam The ParallelTeam instance to execute with.
   */
  public void executeWith(ParallelTeam parallelTeam) {
    try {
      parallelTeam.execute(this);
    } catch (Exception e) {
      String message = " Exception expanding initializing GK.\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  public void init(
      GeneralizedKirkwood generalizedKirkwood,
      Atom[] atoms,
      boolean lambdaTerm,
      AtomicDoubleArray3D grad,
      AtomicDoubleArray3D torque,
      AtomicDoubleArray sharedBornGrad) {
    this.generalizedKirkwood = generalizedKirkwood;
    this.atoms = atoms;
    this.lambdaTerm = lambdaTerm;
    this.grad = grad;
    this.torque = torque;
    this.sharedBornGrad = sharedBornGrad;
  }

  @Override
  public void run() {
    int threadIndex = getThreadIndex();
    if (initializationLoop[threadIndex] == null) {
      initializationLoop[threadIndex] = new InitializationLoop();
    }
    try {
      int nAtoms = atoms.length;
      execute(0, nAtoms - 1, initializationLoop[threadIndex]);
    } catch (Exception e) {
      String message = "Fatal exception initializing coordinates in thread: " + threadIndex + "\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  private class InitializationLoop extends IntegerForLoop {

    private int threadID;

    @Override
    public void run(int lb, int ub) {
      grad.reset(threadID, lb, ub);
      torque.reset(threadID, lb, ub);
      sharedBornGrad.reset(threadID, lb, ub);
      if (lambdaTerm) {
        for (int i = lb; i <= ub; i++) {
          // Update GK parameters.
          generalizedKirkwood.udpateSoluteParameters(i);

          if (!atoms[i].applyLambda()) {
            logger.warning(format(" Atom %s is not alchemical.", atoms[i].toString()));
            logger.warning(" Alchemical GK calculations require all atoms to be alchemical.");
            break;
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
    }
  }
}
