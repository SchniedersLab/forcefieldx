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
package ffx.potential.nonbonded;

import static java.lang.String.format;
import static java.util.Arrays.fill;
import static uk.ac.manchester.tornado.api.collections.math.TornadoMath.abs;
import static uk.ac.manchester.tornado.api.collections.math.TornadoMath.floor;
import static uk.ac.manchester.tornado.api.collections.math.TornadoMath.sqrt;

import ffx.crystal.Crystal;
import ffx.numerics.tornado.FFXTornado;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.VDWType;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import uk.ac.manchester.tornado.api.TaskSchedule;
import uk.ac.manchester.tornado.api.annotations.Parallel;
import uk.ac.manchester.tornado.api.annotations.Reduce;
import uk.ac.manchester.tornado.api.common.TornadoDevice;
import uk.ac.manchester.tornado.api.runtime.TornadoRuntime;

/**
 * The Van der Waals class computes Van der Waals interaction in parallel using a {@link
 * NeighborList} for P1 (no symmetry) {@link Crystal} systems.
 *
 * <p>The repulsive power (e.g. 12), attractive power (e.g. 6) and buffering (e.g. for the AMOEBA
 * buffered-14-7) can all be specified such that both Lennard-Jones and AMOEBA are supported.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class VanDerWaalsTornado extends VanDerWaals {

  private static final Logger logger = Logger.getLogger(VanDerWaalsTornado.class.getName());
  private static final byte XX = 0;
  private static final byte YY = 1;
  private static final byte ZZ = 2;
  /** Lennard-Jones or Buffered 14-7 vdW form. */
  private final VanDerWaalsForm vdwForm;
  /** Cut-off switch. */
  private final double vdwTaper;

  private final double vdwCutoff;
  // *************************************************************************
  // Accumulation variables.
  private int interactions;
  private double energy;
  private double[] grad;
  private Crystal crystal;
  /** An array of all atoms in the system. */
  private Atom[] atoms;
  /** The Force Field that defines the Van der Waals interactions. */
  private ForceField forceField;
  /** A local convenience variable equal to atoms.length. */
  private int nAtoms;
  /** A local reference to the atom class of each atom in the system. */
  private int[] atomClass;
  /** A local copy of atomic coordinates, including reductions on the hydrogen atoms. */
  private double[] coordinates;
  /**
   * Each hydrogen vdW site is located a fraction of the way from the heavy atom nucleus to the
   * hydrogen nucleus (~0.9).
   */
  private double[] reductionValue;
  /**
   * Hydrogen atom vdW sites are located toward their heavy atom relative to their nucleus. This is
   * a look-up that gives the heavy atom index for each hydrogen.
   */
  private int[] reductionIndex;
  /** 1-2, 1-3, 1-4 interactions are masked. */
  private int[] mask;
  /** Pointer into the mask array for each atom. */
  private int[] maskPointer;

  /**
   * The VanDerWaalsTornado class constructor.
   *
   * @param atoms Atom array to do Van Der Waals calculations on.
   * @param crystal The periodic boundary conditions information.
   * @param forceField the ForceField parameters to apply.
   * @param vdwCutoff vdW cutoff.
   * @since 1.0
   */
  public VanDerWaalsTornado(
      Atom[] atoms, Crystal crystal, ForceField forceField, double vdwCutoff) {
    this.atoms = atoms;
    this.crystal = crystal;
    this.forceField = forceField;
    nAtoms = atoms.length;
    vdwForm = new VanDerWaalsForm(forceField);

    // Allocate coordinate arrays and set up reduction indices and values.
    initAtomArrays();

    /*
     Define the multiplicative switch, which sets vdW energy to zero
     at the cutoff distance using a window that begin at 90% of the
     cutoff distance.
    */
    this.vdwCutoff = vdwCutoff;
    this.vdwTaper = 0.9 * vdwCutoff;
    logger.info(toString());
  }

  /** Currently does not use "neighbor-lists" and so is truly N^2. */
  private static void tornadoEnergy(
      int[] atomClass,
      double[] eps,
      double[] rMin,
      double[] reducedXYZ,
      int[] reductionIndex,
      double[] reductionValue,
      double[] bondedScaleFactors,
      int[] maskPointers,
      int[] masks,
      double[] A,
      double[] Ai,
      double[] cutoffs,
      @Reduce double[] energy,
      @Reduce int[] interactions,
      @Reduce double[] grad) {

    // Matrices to apply the minimum image convention.
    // The columns of A are the reciprocal basis vectors
    double A00 = A[0];
    double A01 = A[1];
    double A02 = A[2];
    double A10 = A[3];
    double A11 = A[4];
    double A12 = A[5];
    double A20 = A[6];
    double A21 = A[7];
    double A22 = A[8];
    // a is the first row of A^(-1).
    double Ai00 = Ai[0];
    double Ai01 = Ai[1];
    double Ai02 = Ai[2];
    // b is the second row of A^(-1).
    double Ai10 = Ai[3];
    double Ai11 = Ai[4];
    double Ai12 = Ai[5];
    // c is the third row of A^(-1).
    double Ai20 = Ai[6];
    double Ai21 = Ai[7];
    double Ai22 = Ai[8];

    // Scale factors for bonded atoms.
    double scale12 = bondedScaleFactors[0];
    double scale13 = bondedScaleFactors[1];
    double scale14 = bondedScaleFactors[2];

    // Periodic Boundary Conditions
    boolean aperiodic = false;
    if (cutoffs[0] > 0) {
      aperiodic = true;
    }

    double vdwTaper = cutoffs[1];
    double vdwCutoff = cutoffs[2];
    double vdwTaper2 = vdwTaper * vdwTaper;
    double vdwCutoff2 = vdwCutoff * vdwCutoff;
    boolean gradient = false;
    if (cutoffs[3] > 0) {
      gradient = true;
    }

    // Multiplicative Switch
    double a = vdwTaper;
    double b = vdwCutoff;
    double a2 = a * a;
    double b2 = b * b;
    double ba = b - a;
    double ba2 = ba * ba;
    double denom = ba * ba2 * ba2;
    double c0 = b * b2 * (b2 - 5.0 * a * b + 10.0 * a2) / denom;
    double c1 = -30.0 * a2 * b2 / denom;
    double c2 = 30.0 * b * a * (b + a) / denom;
    double c3 = -10.0 * (a2 + 4.0 * a * b + b2) / denom;
    double c4 = 15.0 * (a + b) / denom;
    double c5 = -6.0 / denom;
    double twoC2 = 2.0 * c2;
    double threeC3 = 3.0 * c3;
    double fourC4 = 4.0 * c4;
    double fiveC5 = 5.0 * c5;

    // Interactions between 1-2, 1-3, 1-4 atoms are masked.
    final int nAtoms = atomClass.length;
    double[] mask = new double[nAtoms];
    for (int i = 0; i < nAtoms; i++) {
      mask[i] = 1.0;
    }

    // AMOEBA Buffered 14-7 parameters.
    final double delta = 0.07;
    final double gamma = 0.12;
    final double delta1 = 1.0 + delta;
    final double d2 = delta1 * delta1;
    final double d4 = d2 * d2;
    final double t1n = delta1 * d2 * d4;
    final double gamma1 = 1.0 + gamma;

    final int XX = 0;
    final int YY = 1;
    final int ZZ = 2;
    // Outer loop over all atoms.
    for (@Parallel int i = 0; i < nAtoms - 1; i++) {
      final int i3 = i * 3;
      final double xi = reducedXYZ[i3 + XX];
      final double yi = reducedXYZ[i3 + YY];
      final double zi = reducedXYZ[i3 + ZZ];
      final int redi = reductionIndex[i];
      final double redv = reductionValue[i];
      final double rediv = 1.0 - redv;
      final int classI = atomClass[i];
      final double ei = eps[classI];
      final double sei = sqrt(ei);
      final double ri = rMin[classI];
      if (ri <= 0.0) {
        continue;
      }
      double gxi = 0.0;
      double gyi = 0.0;
      double gzi = 0.0;
      double gxredi = 0.0;
      double gyredi = 0.0;
      double gzredi = 0.0;

      // Apply masks
      for (int ii = maskPointers[i3]; ii < maskPointers[i3 + 1]; ii++) {
        mask[masks[ii]] = scale12;
      }
      for (int ii = maskPointers[i3 + 1]; ii < maskPointers[i3 + 2]; ii++) {
        mask[masks[ii]] = scale13;
      }
      for (int ii = maskPointers[i3 + 2]; ii < maskPointers[i3 + 3]; ii++) {
        mask[masks[ii]] = scale14;
      }

      // Inner loop over atoms.
      for (int k = i + 1; k < nAtoms; k++) {
        final int k3 = k * 3;
        final double xk = reducedXYZ[k3 + XX];
        final double yk = reducedXYZ[k3 + YY];
        final double zk = reducedXYZ[k3 + ZZ];
        final double[] dx = new double[3];
        dx[0] = xi - xk;
        dx[1] = yi - yk;
        dx[2] = zi - zk;
        double x = dx[0];
        double y = dx[1];
        double z = dx[2];
        double r2;
        if (aperiodic) {
          r2 = x * x + y * y + z * z;
        } else {
          double xf = x * A00 + y * A10 + z * A20;
          double yf = x * A01 + y * A11 + z * A21;
          double zf = x * A02 + y * A12 + z * A22;
          // signum:
          // zero if the argument is zero,
          // 1.0 if the argument is greater than zero,
          // -1.0 if the argument is less than zero.
          // double xfsn = signum(-xf);
          double xfsn = 0.0;
          if (-xf > 0.0) {
            xfsn = 1.0;
          } else if (-xf < 0.0) {
            xfsn = -1.0;
          }
          // double yfsn = signum(-yf);
          double yfsn = 0.0;
          if (-yf > 0.0) {
            yfsn = 1.0;
          } else if (-yf < 0.0) {
            yfsn = -1.0;
          }
          // double zfsn = signum(-zf);
          double zfsn = 0.0;
          if (-zf > 0.0) {
            zfsn = 1.0;
          } else if (-zf < 0.0) {
            zfsn = -1.0;
          }
          xf = floor(abs(xf) + 0.5) * xfsn + xf;
          yf = floor(abs(yf) + 0.5) * yfsn + yf;
          zf = floor(abs(zf) + 0.5) * zfsn + zf;
          x = xf * Ai00 + yf * Ai10 + zf * Ai20;
          y = xf * Ai01 + yf * Ai11 + zf * Ai21;
          z = xf * Ai02 + yf * Ai12 + zf * Ai22;
          dx[0] = x;
          dx[1] = y;
          dx[2] = z;
          r2 = x * x + y * y + z * z;
        }
        final int classK = atomClass[k];
        final double rk = rMin[classK];
        if (r2 <= vdwCutoff2 && mask[k] > 0 && rk > 0) {
          double ri2 = ri * ri;
          double ri3 = ri * ri2;
          double rk2 = rk * rk;
          double rk3 = rk * rk2;
          double irv = 1.0 / (2.0 * (ri3 + rk3) / (ri2 + rk2));
          final double r = sqrt(r2);
          /*
           Calculate Van der Waals interaction energy.
           Notation of Schnieders et al. The structure,
           thermodynamics, and solubility of organic
           crystals from simulation with a polarizable force
           field. J. Chem. Theory Comput. 8, 1721–1736 (2012).
          */
          double ek = eps[classK];
          double sek = sqrt(ek);
          double ev = mask[k] * 4.0 * (ei * ek) / ((sei + sek) * (sei + sek));
          final double rho = r * irv;
          final double rho2 = rho * rho;
          final double rhoDisp1 = rho2 * rho2 * rho2;
          final double rhoDisp = rhoDisp1 * rho;
          final double rhoD = rho + delta;
          final double rhoD2 = rhoD * rhoD;
          final double rhoDelta1 = rhoD2 * rhoD2 * rhoD2;
          final double rhoDelta = rhoDelta1 * (rho + delta);
          final double rhoDispGamma = rhoDisp + gamma;
          final double t1d = 1.0 / rhoDelta;
          final double t2d = 1.0 / rhoDispGamma;
          final double t1 = t1n * t1d;
          final double t2a = gamma1 * t2d;
          final double t2 = t2a - 2.0;
          double eik = ev * t1 * t2;
          /*
           Apply a multiplicative switch if the interaction
           distance is greater than the beginning of the taper.
          */
          double taper = 1.0;
          double dtaper = 0.0;
          if (r2 > vdwTaper2) {
            final double r3 = r2 * r;
            final double r4 = r2 * r2;
            final double r5 = r2 * r3;
            taper = c5 * r5 + c4 * r4 + c3 * r3 + c2 * r2 + c1 * r + c0;
            dtaper = fiveC5 * r4 + fourC4 * r3 + threeC3 * r2 + twoC2 * r + c1;
          }
          eik *= taper;
          energy[0] += eik;
          interactions[0] += 1;
          if (!gradient) {
            continue;
          }
          final int redk = reductionIndex[k];
          final double red = reductionValue[k];
          final double redkv = 1.0 - red;
          final double dt1d_dr = 7.0 * rhoDelta1 * irv;
          final double dt2d_dr = 7.0 * rhoDisp1 * irv;
          final double dt1_dr = t1 * dt1d_dr * t1d;
          final double dt2_dr = t2a * dt2d_dr * t2d;
          final double dedr = -ev * (dt1_dr * t2 + t1 * dt2_dr);
          final double ir = 1.0 / r;
          final double drdx = dx[0] * ir;
          final double drdy = dx[1] * ir;
          final double drdz = dx[2] * ir;
          final double dswitch = (eik * dtaper + dedr * taper);
          final double dedx = dswitch * drdx;
          final double dedy = dswitch * drdy;
          final double dedz = dswitch * drdz;
          gxi += dedx * redv;
          gyi += dedy * redv;
          gzi += dedz * redv;
          gxredi += dedx * rediv;
          gyredi += dedy * rediv;
          gzredi += dedz * rediv;
          // Atom K
          grad[k3 + XX] -= red * dedx;
          grad[k3 + YY] -= red * dedy;
          grad[k3 + ZZ] -= red * dedz;
          // Atom K is reduced by Atom redK;
          int r3 = redk * 3;
          grad[r3 + XX] -= redkv * dedx;
          grad[r3 + YY] -= redkv * dedy;
          grad[r3 + ZZ] -= redkv * dedz;
        }
      }
      if (gradient) {
        // Atom I gradient.
        grad[i3 + XX] += gxi;
        grad[i3 + YY] += gyi;
        grad[i3 + ZZ] += gzi;
        // Atom I is reduced by Atom redI;
        int r3 = redi * 3;
        grad[r3 + XX] += gxredi;
        grad[r3 + YY] += gyredi;
        grad[r3 + ZZ] += gzredi;
      }

      // Remove masks
      for (int ii = maskPointers[i3]; ii < maskPointers[i3 + 1]; ii++) {
        mask[masks[ii]] = 1.0;
      }
      for (int ii = maskPointers[i3 + 1]; ii < maskPointers[i3 + 2]; ii++) {
        mask[masks[ii]] = 1.0;
      }
      for (int ii = maskPointers[i3 + 2]; ii < maskPointers[i3 + 3]; ii++) {
        mask[masks[ii]] = 1.0;
      }
    }
  }

  /**
   * The energy routine may be called repeatedly.
   *
   * @param gradient If true, gradients with respect to atomic coordinates are computed.
   * @param print If true, there is verbose printing.
   * @return The energy.
   * @since 1.0
   */
  public double energy(boolean gradient, boolean print) {

    if (vdwForm.vdwType != VanDerWaalsForm.VDW_TYPE.BUFFERED_14_7) {
      logger.severe((" TornadoVM vdW only supports AMOEBA."));
    }

    // Initialize coordinates and gradient array.
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      // Load atomic coordinates.
      double x = atom.getX();
      double y = atom.getY();
      double z = atom.getZ();
      int i3 = i * 3;
      coordinates[i3 + XX] = x;
      coordinates[i3 + YY] = y;
      coordinates[i3 + ZZ] = z;
    }

    double[] eps = vdwForm.getEps();
    double[] rmin = vdwForm.getRmin();

    // Compute reduced hydrogen atom coordinates.
    double[] reducedXYZ = new double[nAtoms * 3];
    final byte XX = 0;
    final byte YY = 1;
    final byte ZZ = 2;
    for (int i = 0; i < nAtoms; i++) {
      int i3 = i * 3;
      double x = coordinates[i3 + XX];
      double y = coordinates[i3 + YY];
      double z = coordinates[i3 + ZZ];
      int redIndex = reductionIndex[i];
      if (redIndex >= 0) {
        int r3 = redIndex * 3;
        double rx = coordinates[r3 + XX];
        double ry = coordinates[r3 + YY];
        double rz = coordinates[r3 + ZZ];
        double r = reductionValue[i];
        reducedXYZ[i3 + XX] = r * (x - rx) + rx;
        reducedXYZ[i3 + YY] = r * (y - ry) + ry;
        reducedXYZ[i3 + ZZ] = r * (z - rz) + rz;
      } else {
        reducedXYZ[i3 + XX] = x;
        reducedXYZ[i3 + YY] = y;
        reducedXYZ[i3 + ZZ] = z;
      }
    }

    // Gradient flag.
    double doGradient = 0.0;
    if (gradient) {
      doGradient = 1.0;
    }

    // Initialize periodic boundary conditions.
    Crystal c = crystal;
    double[] A = {c.A00, c.A01, c.A02, c.A10, c.A11, c.A12, c.A20, c.A21, c.A22};
    double[] Ai = {c.Ai00, c.Ai01, c.Ai02, c.Ai10, c.Ai11, c.Ai12, c.Ai20, c.Ai21, c.Ai22};
    double[] bondedScaleFactors = {vdwForm.scale12, vdwForm.scale13, vdwForm.scale14};
    double aperiodic = 0.0;
    if (crystal.aperiodic()) {
      aperiodic = 1.0;
    }

    // Periodic boundary information, plus include the gradient flag.
    double[] cutoffs = {aperiodic, vdwTaper, vdwCutoff, doGradient};

    // Output
    double[] energy = new double[1];
    int[] interactions = new int[1];
    if (gradient) {
      fill(grad, 0.0);
    }

    // Check the method without TornadoVM
    tornadoEnergy(
        atomClass,
        eps,
        rmin,
        reducedXYZ,
        reductionIndex,
        reductionValue,
        bondedScaleFactors,
        maskPointer,
        mask,
        A,
        Ai,
        cutoffs,
        energy,
        interactions,
        grad);

    logger.info(format(" JVM: %16.8f %d", energy[0], interactions[0]));

    energy[0] = 0.0;
    interactions[0] = 0;
    if (gradient) {
      fill(grad, 0.0);
    }

    TornadoDevice device = TornadoRuntime.getTornadoRuntime().getDefaultDevice();
    FFXTornado.logDevice(device);
    TaskSchedule graph =
        new TaskSchedule("vdW")
            .streamIn(
                atomClass,
                eps,
                rmin,
                reducedXYZ,
                reductionIndex,
                reductionValue,
                bondedScaleFactors,
                maskPointer,
                mask,
                A,
                Ai,
                cutoffs,
                energy,
                interactions,
                grad)
            .task(
                "energy",
                VanDerWaalsTornado::tornadoEnergy,
                atomClass,
                eps,
                rmin,
                reducedXYZ,
                reductionIndex,
                reductionValue,
                bondedScaleFactors,
                maskPointer,
                mask,
                A,
                Ai,
                cutoffs,
                energy,
                interactions,
                grad)
            .streamOut(energy, interactions, grad);

    graph.setDevice(device);
    graph.warmup();
    graph.execute();
    graph.dumpProfiles();
    device.reset();

    logger.info(format(" Tornado OpenCL: %16.8f %d", energy[0], interactions[0]));

    // Add gradient to atoms.
    if (gradient) {
      for (int i = 0; i < nAtoms - 1; i++) {
        Atom ai = atoms[i];
        int i3 = i * 3;
        ai.addToXYZGradient(grad[i3 + XX], grad[i3 + YY], grad[i3 + ZZ]);
      }
    }

    this.energy = energy[0];
    this.interactions = interactions[0];
    return this.energy;
  }

  /**
   * Get the total Van der Waals potential energy.
   *
   * @return The energy.
   * @since 1.0
   */
  public double getEnergy() {
    return energy;
  }

  /**
   * Get the number of interacting pairs.
   *
   * @return The interaction count.
   * @since 1.0
   */
  public int getInteractions() {
    return interactions;
  }

  /**
   * Setter for the field <code>atoms</code>.
   *
   * @param atoms an array of {@link Atom} objects.
   */
  public void setAtoms(Atom[] atoms) {
    this.atoms = atoms;
    this.nAtoms = atoms.length;
    initAtomArrays();
  }

  /**
   * If the crystal being passed in is not equal to the current crystal, then some Van der Waals
   * data structures may need to updated. If <code>nSymm</code> has changed, update arrays
   * dimensioned by nSymm. Finally, rebuild the neighbor-lists.
   *
   * @param crystal The new crystal instance defining the symmetry and boundary conditions.
   */
  public void setCrystal(Crystal crystal) {
    this.crystal = crystal;
    int newNSymm = crystal.getNumSymOps();
    if (newNSymm != 1) {
      String message = " SymOps are not supported by VanDerWaalsTornado.\n";
      logger.log(Level.SEVERE, message);
    }
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    StringBuffer sb = new StringBuffer("\n  Van der Waals\n");
    sb.append(format("   Switch Start:                         %6.3f (A)\n", vdwTaper));
    sb.append(format("   Cut-Off:                              %6.3f (A)\n", vdwCutoff));
    return sb.toString();
  }

  /** Allocate coordinate arrays and set up reduction indices and values. */
  private void initAtomArrays() {
    if (atomClass == null || nAtoms > atomClass.length) {
      atomClass = new int[nAtoms];
      coordinates = new double[nAtoms * 3];
      reductionIndex = new int[nAtoms];
      reductionValue = new double[nAtoms];
      grad = new double[nAtoms * 3];
      maskPointer = new int[nAtoms * 3 + 1];
    }

    int numBonds = 0;
    int numAngles = 0;
    int numTorsions = 0;
    for (int i = 0; i < nAtoms; i++) {
      Atom ai = atoms[i];
      numBonds += ai.getNumBonds();
      numAngles += ai.getNumAngles();
      numTorsions += ai.getNumDihedrals();
    }
    mask = new int[numBonds + numAngles + numTorsions];

    int[][] mask12 = getMask12();
    int[][] mask13 = getMask13();
    int[][] mask14 = getMask14();

    int index = 0;
    for (int i = 0; i < nAtoms; i++) {
      Atom ai = atoms[i];
      assert (i == ai.getXyzIndex() - 1);
      double[] xyz = ai.getXYZ(null);
      int i3 = i * 3;
      coordinates[i3 + XX] = xyz[XX];
      coordinates[i3 + YY] = xyz[YY];
      coordinates[i3 + ZZ] = xyz[ZZ];
      AtomType atomType = ai.getAtomType();
      if (atomType == null) {
        logger.severe(ai.toString());
      }
      String vdwIndex = forceField.getString("VDWINDEX", "Class");
      if (vdwIndex.equalsIgnoreCase("Type")) {
        atomClass[i] = atomType.type;
      } else {
        atomClass[i] = atomType.atomClass;
      }
      VDWType type = forceField.getVDWType(Integer.toString(atomClass[i]));
      if (type == null) {
        logger.info(" No VdW type for atom class " + atomClass[i]);
        logger.severe(" No VdW type for atom " + ai.toString());
        return;
      }
      ai.setVDWType(type);
      List<Bond> bonds = ai.getBonds();
      numBonds = bonds.size();
      if (type.reductionFactor > 0.0 && numBonds == 1) {
        Bond bond = bonds.get(0);
        Atom heavyAtom = bond.get1_2(ai);
        // Atom indexes start at 1
        reductionIndex[i] = heavyAtom.getIndex() - 1;
        reductionValue[i] = type.reductionFactor;
      } else {
        reductionIndex[i] = i;
        reductionValue[i] = 0.0;
      }
      // Store bond mask for atom i.
      maskPointer[3 * i] = index;
      for (int value : mask12[i]) {
        mask[index++] = value;
      }
      maskPointer[3 * i + 1] = index;
      for (int value : mask13[i]) {
        mask[index++] = value;
      }
      maskPointer[3 * i + 2] = index;
      for (int value : mask14[i]) {
        mask[index++] = value;
      }
    }
    maskPointer[3 * nAtoms] = index;
  }

  /**
   * Log the Van der Waals interaction.
   *
   * @param i Atom i.
   * @param k Atom j.
   * @param minr The minimum vdW separation distance.
   * @param r The distance rij.
   * @param eij The interaction energy.
   * @since 1.0
   */
  private void log(int i, int k, double minr, double r, double eij) {
    logger.info(
        format(
            "VDW %6d-%s %6d-%s %10.4f  %10.4f  %10.4f",
            atoms[i].getIndex(),
            atoms[i].getAtomType().name,
            atoms[k].getIndex(),
            atoms[k].getAtomType().name,
            1.0 / minr,
            r,
            eij));
  }
}
