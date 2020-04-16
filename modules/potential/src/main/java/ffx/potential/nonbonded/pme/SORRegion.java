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

import static ffx.potential.parameters.MultipoleType.t001;
import static ffx.potential.parameters.MultipoleType.t010;
import static ffx.potential.parameters.MultipoleType.t100;
import static org.apache.commons.math3.util.FastMath.max;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.reduction.SharedDouble;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.ParticleMeshEwaldCart.EwaldParameters;
import ffx.potential.parameters.ForceField;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Parallel successive over-relaxation (SOR) solver for the self-consistent field.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class SORRegion extends ParallelRegion {

  private static final Logger logger = Logger.getLogger(SORRegion.class.getName());
  private final int maxThreads;
  private final double polsor;
  private final SORLoop[] sorLoop;
  private final SharedDouble sharedEps;
  private final SharedDouble sharedEpsCR;
  /** Dimensions of [nsymm][nAtoms][3] */
  public double[][][] inducedDipole;

  public double[][][] inducedDipoleCR;
  /** Direct induced dipoles. */
  public double[][] directDipole;

  public double[][] directDipoleCR;
  /** An ordered array of atoms in the system. */
  private Atom[] atoms;

  private double[] polarizability;
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

  public SORRegion(int nt, ForceField forceField) {
    maxThreads = nt;
    sorLoop = new SORLoop[nt];
    sharedEps = new SharedDouble();
    sharedEpsCR = new SharedDouble();
    polsor = forceField.getDouble("POLAR_SOR", 0.70);
  }

  public double getEps() {
    double eps = sharedEps.get();
    double epsCR = sharedEpsCR.get();
    return max(eps, epsCR);
  }

  public double getSOR() {
    return polsor;
  }

  public void init(
      Atom[] atoms,
      double[] polarizability,
      double[][][] inducedDipole,
      double[][][] inducedDipoleCR,
      double[][] directDipole,
      double[][] directDipoleCR,
      double[][] cartesianDipolePhi,
      double[][] cartesianDipolePhiCR,
      AtomicDoubleArray3D field,
      AtomicDoubleArray3D fieldCR,
      boolean generalizedKirkwoodTerm,
      GeneralizedKirkwood generalizedKirkwood,
      EwaldParameters ewaldParameters) {
    this.atoms = atoms;
    this.polarizability = polarizability;
    this.inducedDipole = inducedDipole;
    this.inducedDipoleCR = inducedDipoleCR;
    this.directDipole = directDipole;
    this.directDipoleCR = directDipoleCR;
    this.cartesianDipolePhi = cartesianDipolePhi;
    this.cartesianDipolePhiCR = cartesianDipolePhiCR;
    this.field = field;
    this.fieldCR = fieldCR;
    this.generalizedKirkwoodTerm = generalizedKirkwoodTerm;
    this.generalizedKirkwood = generalizedKirkwood;
    this.aewald = ewaldParameters.aewald;
    this.aewald3 = ewaldParameters.aewald3;
  }

  @Override
  public void run() throws Exception {
    try {
      int ti = getThreadIndex();
      if (sorLoop[ti] == null) {
        sorLoop[ti] = new SORLoop();
      }
      int nAtoms = atoms.length;
      execute(0, nAtoms - 1, sorLoop[ti]);
    } catch (RuntimeException ex) {
      logger.warning(
          "Fatal exception computing the mutual induced dipoles in thread " + getThreadIndex());
      throw ex;
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
    sharedEps.set(0.0);
    sharedEpsCR.set(0.0);
  }

  private class SORLoop extends IntegerForLoop {

    private double eps, epsCR;

    @Override
    public void finish() {
      sharedEps.addAndGet(eps);
      sharedEpsCR.addAndGet(epsCR);
    }

    @Override
    public void run(int lb, int ub) throws Exception {
      int threadID = getThreadIndex();
      final double[][] induced0 = inducedDipole[0];
      final double[][] inducedCR0 = inducedDipoleCR[0];
      if (aewald > 0.0) {
        // Add the self and reciprocal space fields to the real space field.
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

      if (generalizedKirkwoodTerm) {
        AtomicDoubleArray3D fieldGK = generalizedKirkwood.getFieldGK();
        AtomicDoubleArray3D fieldGKCR = generalizedKirkwood.getFieldGKCR();
        // Add the GK reaction field to the intra-molecular field.
        for (int i = lb; i <= ub; i++) {
          field.add(threadID, i, fieldGK.getX(i), fieldGK.getY(i), fieldGK.getZ(i));
          fieldCR.add(threadID, i, fieldGKCR.getX(i), fieldGKCR.getY(i), fieldGKCR.getZ(i));
        }
      }

      // Reduce the real space field.
      field.reduce(lb, ub);
      fieldCR.reduce(lb, ub);

      // Apply Successive Over-Relaxation (SOR).
      for (int i = lb; i <= ub; i++) {
        final double[] ind = induced0[i];
        final double[] indCR = inducedCR0[i];
        final double[] direct = directDipole[i];
        final double[] directCR = directDipoleCR[i];
        final double polar = polarizability[i];
        for (int j = 0; j < 3; j++) {
          double previous = ind[j];
          double mutual = polar * field.get(j, i);
          ind[j] = direct[j] + mutual;
          double delta = polsor * (ind[j] - previous);
          ind[j] = previous + delta;
          eps += delta * delta;
          previous = indCR[j];
          mutual = polar * fieldCR.get(j, i);
          indCR[j] = directCR[j] + mutual;
          delta = polsor * (indCR[j] - previous);
          indCR[j] = previous + delta;
          epsCR += delta * delta;
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
