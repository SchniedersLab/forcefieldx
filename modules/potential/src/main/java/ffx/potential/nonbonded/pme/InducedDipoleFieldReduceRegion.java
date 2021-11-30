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

import static ffx.potential.parameters.MultipoleType.t001;
import static ffx.potential.parameters.MultipoleType.t010;
import static ffx.potential.parameters.MultipoleType.t100;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.ParticleMeshEwaldCart.EwaldParameters;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Parallel summation and reduction of components of the induced dipole field at each atom.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class InducedDipoleFieldReduceRegion extends ParallelRegion {

  private static final Logger logger = Logger.getLogger(DirectRegion.class.getName());
  private final InducedDipoleFieldReduceLoop[] inducedDipoleFieldReduceLoop;
  /** Dimensions of [nsymm][nAtoms][3] */
  public double[][][] inducedDipole;

  public double[][][] inducedDipoleCR;
  /** An ordered array of atoms in the system. */
  private Atom[] atoms;

  private double[][] cartesianDipolePhi;
  private double[][] cartesianDipolePhiCR;
  /** Field array. */
  private AtomicDoubleArray3D field;
  /** Chain rule field array. */
  private AtomicDoubleArray3D fieldCR;
  /** Flag to indicate use of generalized Kirkwood. */
  private boolean generalizedKirkwoodTerm;

  private GeneralizedKirkwood generalizedKirkwood;
  private double aewald;
  private double aewald3;

  public InducedDipoleFieldReduceRegion(int nt) {
    inducedDipoleFieldReduceLoop = new InducedDipoleFieldReduceLoop[nt];
  }

  /**
   * Execute the InducedDipoleFieldReduceRegion with the passed ParallelTeam.
   *
   * @param parallelTeam The ParallelTeam instance to execute with.
   */
  public void executeWith(ParallelTeam parallelTeam) {
    try {
      parallelTeam.execute(this);
    } catch (Exception e) {
      String message = " Exception computing induced dipole field.\n";
      logger.log(Level.WARNING, message, e);
    }
  }

  public void init(
      Atom[] atoms,
      double[][][] inducedDipole,
      double[][][] inducedDipoleCR,
      boolean generalizedKirkwoodTerm,
      GeneralizedKirkwood generalizedKirkwood,
      EwaldParameters ewaldParameters,
      double[][] cartesianDipolePhi,
      double[][] cartesianDipolePhiCR,
      AtomicDoubleArray3D field,
      AtomicDoubleArray3D fieldCR) {
    // Input
    this.atoms = atoms;
    this.inducedDipole = inducedDipole;
    this.inducedDipoleCR = inducedDipoleCR;
    this.generalizedKirkwoodTerm = generalizedKirkwoodTerm;
    this.generalizedKirkwood = generalizedKirkwood;
    this.aewald = ewaldParameters.aewald;
    this.aewald3 = ewaldParameters.aewald3;
    this.cartesianDipolePhi = cartesianDipolePhi;
    this.cartesianDipolePhiCR = cartesianDipolePhiCR;
    // Output
    this.field = field;
    this.fieldCR = fieldCR;
  }

  @Override
  public void run() throws Exception {
    try {
      int threadID = getThreadIndex();
      if (inducedDipoleFieldReduceLoop[threadID] == null) {
        inducedDipoleFieldReduceLoop[threadID] = new InducedDipoleFieldReduceLoop();
      }
      int nAtoms = atoms.length;
      execute(0, nAtoms - 1, inducedDipoleFieldReduceLoop[threadID]);
    } catch (Exception e) {
      String message =
          " Fatal exception computing the mutual induced dipoles in thread "
              + getThreadIndex()
              + "\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  private class InducedDipoleFieldReduceLoop extends IntegerForLoop {

    @Override
    public void run(int lb, int ub) throws Exception {
      int threadID = getThreadIndex();
      final double[][] induced0 = inducedDipole[0];
      final double[][] inducedCR0 = inducedDipoleCR[0];
      // Add the PME self and reciprocal space fields to the real space field.
      if (aewald > 0.0) {
        for (int i = lb; i <= ub; i++) {
          double[] dipolei = induced0[i];
          double[] dipoleCRi = inducedCR0[i];
          final double[] phii = cartesianDipolePhi[i];
          final double[] phiCRi = cartesianDipolePhiCR[i];
          double fx = aewald3 * dipolei[0] - phii[t100];
          double fy = aewald3 * dipolei[1] - phii[t010];
          double fz = aewald3 * dipolei[2] - phii[t001];
          double fxCR = aewald3 * dipoleCRi[0] - phiCRi[t100];
          double fyCR = aewald3 * dipoleCRi[1] - phiCRi[t010];
          double fzCR = aewald3 * dipoleCRi[2] - phiCRi[t001];
          field.add(threadID, i, fx, fy, fz);
          fieldCR.add(threadID, i, fxCR, fyCR, fzCR);
        }
      }
      // Add the GK reaction field.
      if (generalizedKirkwoodTerm) {
        AtomicDoubleArray3D fieldGK = generalizedKirkwood.getFieldGK();
        AtomicDoubleArray3D fieldGKCR = generalizedKirkwood.getFieldGKCR();
        for (int i = lb; i <= ub; i++) {
          field.add(threadID, i, fieldGK.getX(i), fieldGK.getY(i), fieldGK.getZ(i));
          fieldCR.add(threadID, i, fieldGKCR.getX(i), fieldGKCR.getY(i), fieldGKCR.getZ(i));
        }
      }
      // Reduce the PME and GK contributions.
      field.reduce(lb, ub);
      fieldCR.reduce(lb, ub);
    }

    @Override
    public IntegerSchedule schedule() {
      return IntegerSchedule.fixed();
    }
  }
}
