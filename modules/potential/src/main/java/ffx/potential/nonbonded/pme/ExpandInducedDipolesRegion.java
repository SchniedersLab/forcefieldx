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
package ffx.potential.nonbonded.pme;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Atom;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Parallel expansion of the asymmetric unit induced dipoles to symmetry mates by applying symmetry
 * operator rotation matrices.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ExpandInducedDipolesRegion extends ParallelRegion {

  private static final Logger logger = Logger.getLogger(ExpandInducedDipolesRegion.class.getName());
  private final ExpandInducedDipoleLoop[] expandInducedDipoleLoop;
  /** Dimensions of [nsymm][nAtoms][3] */
  public double[][][] inducedDipole;

  public double[][][] inducedDipoleCR;
  /** An ordered array of atoms in the system. */
  private Atom[] atoms;
  /** Unit cell and spacegroup information. */
  private Crystal crystal;

  public ExpandInducedDipolesRegion(int maxThreads) {
    expandInducedDipoleLoop = new ExpandInducedDipoleLoop[maxThreads];
    for (int i = 0; i < maxThreads; i++) {
      expandInducedDipoleLoop[i] = new ExpandInducedDipoleLoop();
    }
  }

  /**
   * Execute the ExpandInducedDipolesRegion with the passed ParallelTeam.
   *
   * @param parallelTeam The ParallelTeam instance to execute with.
   */
  public void executeWith(ParallelTeam parallelTeam) {
    try {
      parallelTeam.execute(this);
    } catch (Exception e) {
      String message = " Exception expanding induced dipoles.\n";
      logger.log(Level.WARNING, message, e);
    }
  }

  public void init(
      Atom[] atoms, Crystal crystal, double[][][] inducedDipole, double[][][] inducedDipoleCR) {
    // Input
    this.atoms = atoms;
    this.crystal = crystal;
    // Output
    this.inducedDipole = inducedDipole;
    this.inducedDipoleCR = inducedDipoleCR;
  }

  @Override
  public void run() {
    try {
      int nAtoms = atoms.length;
      execute(0, nAtoms - 1, expandInducedDipoleLoop[getThreadIndex()]);
    } catch (Exception e) {
      String message =
          "Fatal exception expanding coordinates in thread: " + getThreadIndex() + "\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  private class ExpandInducedDipoleLoop extends IntegerForLoop {

    @Override
    public void run(int lb, int ub) {
      List<SymOp> symOps = crystal.spaceGroup.symOps;
      int nSymm = symOps.size();
      for (int s = 1; s < nSymm; s++) {
        SymOp symOp = crystal.spaceGroup.symOps.get(s);
        for (int ii = lb; ii <= ub; ii++) {
          crystal.applySymRot(inducedDipole[0][ii], inducedDipole[s][ii], symOp);
          crystal.applySymRot(inducedDipoleCR[0][ii], inducedDipoleCR[s][ii], symOp);
        }
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return IntegerSchedule.fixed();
    }
  }
}
