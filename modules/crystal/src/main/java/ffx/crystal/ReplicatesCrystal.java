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
package ffx.crystal;

import static java.lang.String.format;

import java.util.List;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The ReplicatesCrystal class extends Crystal to generate additional symmetry operators needed to
 * describe a "replicated" super cell. The replicated crystal cell edges are of length {l*a, m*b,
 * n*c} where l, m and n are integers and a, b and c are the original unit cell edge lengths. <br>
 * The replicates integers l, m and n are chosen large enough for the ReplicatesCrystal to allow
 * consistent application of the minimum image convention. This is ensured by increasing l, m and/or
 * n until a sphere of of necessary radius fits entirely inside the ReplicatedCrystal. <br>
 *
 * @author Michael J. Schnieders
 * @see Crystal
 * @since 1.0
 */
public class ReplicatesCrystal extends Crystal {

  /** The logger. */
  private static final Logger logger = Logger.getLogger(ReplicatesCrystal.class.getName());
  /** The base unit cell for the system being simulated. */
  private final Crystal unitCell;
  /** The cut-off distance in Angstroms. */
  private final double cutOff;
  /** The number of replicates along the a-axis. */
  private int l;
  /** The number of replicates along the b-axis. */
  private int m;
  /** The number of replicates along the c-axis. */
  private int n;

  /**
   * Constructor for a ReplicatesCrystal.
   *
   * @param unitCell The base unit cell.
   * @param l Number of replicates along the a-axis.
   * @param m Number of replicates along the b-axis.
   * @param n Number of replicates along the c-axis.
   * @param cutOff2 Twice the cut-off distance.
   * @since 1.0
   */
  public ReplicatesCrystal(Crystal unitCell, int l, int m, int n, double cutOff2) {
    super(
        unitCell.a * l,
        unitCell.b * m,
        unitCell.c * n,
        unitCell.alpha,
        unitCell.beta,
        unitCell.gamma,
        unitCell.spaceGroup.shortName);
    this.unitCell = unitCell;

    assert (l >= 1);
    assert (m >= 1);
    assert (n >= 1);
    this.l = l;
    this.m = m;
    this.n = n;
    this.cutOff = cutOff2 / 2.0;

    /*
     At this point, the ReplicatesCrystal references a SpaceGroup instance
     whose symmetry operators are inconsistent. This is corrected by
     generating symmetry operators to fill up the ReplicatesCrystal based
     on the asymmetric unit.
    */
    updateReplicateOperators();
  }

  /**
   * Returns a ReplicatesCrystal large enough to satisfy the minimum image convention for the
   * specified unit cell and cutoff criteria. If the unit cell is already sufficiently large, then it
   * is returned.
   *
   * @param unitCell The unit cell of the crystal.
   * @param cutOff2 Two times the cutoff distance.
   * @return A Crystal or ReplicatesCrystal large enough to satisfy the minimum image convention.
   */
  public static Crystal replicatesCrystalFactory(Crystal unitCell, double cutOff2) {

    if (unitCell == null || unitCell.aperiodic()) {
      return unitCell;
    }

    int l = 1;
    int m = 1;
    int n = 1;

    double cutOff = cutOff2 / 2.0;

    while (unitCell.interfacialRadiusA * l < cutOff) {
      l++;
    }
    while (unitCell.interfacialRadiusB * m < cutOff) {
      m++;
    }
    while (unitCell.interfacialRadiusC * n < cutOff) {
      n++;
    }

    if (l * m * n > 1) {
      return new ReplicatesCrystal(unitCell, l, m, n, cutOff2);
    } else {
      return unitCell;
    }
  }

  /**
   * Change the cell vectors for the base unit cell, which is followed by an update of the
   * ReplicateCrystal parameters and possibly the number of replicated cells.
   *
   * @param cellVectors 3x3 matrix of cell vectors.
   * @return True if the perturbation of cell vectors succeeds.
   */
  public boolean setCellVectors(double[][] cellVectors) {

    // First, update the parameters of the unit cell.
    if (unitCell.setCellVectors(cellVectors)) {

      // Then, update the parameters of the ReplicatesCrystal and possibly the number of replicates.
      return updateReplicatesDimensions();
    }
    return false;
  }

  /**
   * Change the cell vectors and volume for the base unit cell, which is followed by an update of the
   * ReplicateCrystal parameters and possibly the number of replicated cells.
   *
   * @param cellVectors 3x3 matrix of cell vectors.
   * @param targetAUVolume the target volume for the new cell Vectors.
   * @return True if the perturbation of cell vectors succeeds.
   */
  public boolean setCellVectorsAndVolume(double[][] cellVectors, double targetAUVolume) {
    // First, update the parameters of the unit cell.
    if (unitCell.setCellVectorsAndVolume(cellVectors, targetAUVolume)) {

      // Then, update the parameters of the ReplicatesCrystal and possibly the number of replicates.
      return updateReplicatesDimensions();
    }
    return false;
  }

  /**
   * Change the cell parameters for the base unit cell, which is followed by an update of the
   * ReplicateCrystal parameters and possibly the number of replicated cells.
   *
   * @param a The length of the a-axis for the base unit cell (in Angstroms).
   * @param b The length of the b-axis for the base unit cell (in Angstroms).
   * @param c The length of the c-axis for the base unit cell (in Angstroms).
   * @param alpha The angle between the b-axis and c-axis (in Degrees).
   * @param beta The angle between the a-axis and c-axis (in Degrees).
   * @param gamma The angle between the a-axis and b-axis (in Degrees).
   * @return True is returned if the unit cell and replicates cell are updated successfully.
   */
  @Override
  public boolean changeUnitCellParameters(
      double a, double b, double c, double alpha, double beta, double gamma) {

    // First, update the parameters of the unit cell.
    if (unitCell.changeUnitCellParameters(a, b, c, alpha, beta, gamma)) {

      // Then, update the parameters of the ReplicatesCrystal and possibly the number of replicates.
      return updateReplicatesDimensions();
    }
    return false;
  }

  /**
   * Change the cell parameters for the base unit cell, which is followed by an update of the
   * ReplicateCrystal parameters and possibly the number of replicated cells.
   *
   * @param a The length of the a-axis for the base unit cell (in Angstroms).
   * @param b The length of the b-axis for the base unit cell (in Angstroms).
   * @param c The length of the c-axis for the base unit cell (in Angstroms).
   * @param alpha The angle between the b-axis and c-axis (in Degrees).
   * @param beta The angle between the a-axis and c-axis (in Degrees).
   * @param gamma The angle between the a-axis and b-axis (in Degrees).
   * @param targetAUVolume the target volume for the new cell Vectors.
   * @return True is returned if the unit cell and replicates cell are updated successfully.
   */
  @Override
  public boolean changeUnitCellParametersAndVolume(
      double a, double b, double c, double alpha, double beta, double gamma, double targetAUVolume) {

    // First, update the parameters of the unit cell.
    if (unitCell.changeUnitCellParametersAndVolume(a, b, c, alpha, beta, gamma, targetAUVolume)) {

      // Then, update the parameters of the ReplicatesCrystal and possibly the number of replicates.
      return updateReplicatesDimensions();
    }
    return false;
  }


  /** Two crystals are equal only if all unit cell parameters are exactly the same. */
  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    ReplicatesCrystal replicatesCrystal = (ReplicatesCrystal) o;
    return (Objects.equals(unitCell, replicatesCrystal.unitCell)
        && a == replicatesCrystal.a
        && b == replicatesCrystal.b
        && c == replicatesCrystal.c
        && alpha == replicatesCrystal.alpha
        && beta == replicatesCrystal.beta
        && gamma == replicatesCrystal.gamma
        && spaceGroup.number == replicatesCrystal.spaceGroup.number
        && spaceGroup.symOps.size() == replicatesCrystal.spaceGroup.symOps.size());
  }

  /** {@inheritDoc} */
  @Override
  public boolean getCheckRestrictions() {
    return unitCell.getCheckRestrictions();
  }

  /** {@inheritDoc} */
  @Override
  public void setCheckRestrictions(boolean checkRestrictions) {
    this.checkRestrictions = checkRestrictions;
    unitCell.setCheckRestrictions(checkRestrictions);
  }

  /**
   * Return the density of the ReplicatesCrystal.
   *
   * @param mass The mass of the ReplicatesCrystal.
   * @return The density.
   */
  public double getDensity(double mass) {
    return unitCell.getDensity(mass);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Returns the unit cell for this ReplicatesCrystal. This is useful for the reciprocal space
   * portion of PME that operates on the unit cell even though the real space cutoff requires a
   * ReplicatesCrystal.
   */
  @Override
  public Crystal getUnitCell() {
    return unitCell;
  }

  /**
   * Update the ReplicatesCrystal using random parameters with the target density.
   *
   * @param density Target density.
   * @param mass Mass of the ReplicatesCrystal.
   */
  public boolean randomParameters(double density, double mass) {
    boolean succeed = unitCell.randomParameters(density, mass);
    if (succeed) {
      updateReplicatesDimensions();
    }
    return succeed;
  }

  /**
   * Update the ReplicatesCrystal dimensions to the target density.
   *
   * @param density Target density.
   * @param mass Mass of the ReplicatesCrystal.
   */
  public void setDensity(double density, double mass) {
    unitCell.setDensity(density, mass);
    updateReplicatesDimensions();
  }

  /**
   * {@inheritDoc}
   *
   * <p>Include information about the base unit cell and replicates cell.
   */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder(unitCell.toString());

    sb.append("\n\n Replicates Cell\n");
    sb.append(format("  Dimension:                    (%3d x%3d x%3d)\n", l, m, n));
    sb.append(format("  A-axis:                              %8.3f\n", a));
    sb.append(format("  B-axis:                              %8.3f\n", b));
    sb.append(format("  C-axis:                              %8.3f\n", c));
    sb.append(format("  Alpha:                               %8.3f\n", alpha));
    sb.append(format("  Beta:                                %8.3f\n", beta));
    sb.append(format("  Gamma:                               %8.3f\n", gamma));
    sb.append(format("  Total Symmetry Operators:            %8d", spaceGroup.getNumberOfSymOps()));
    return sb.toString();
  }

  private boolean updateReplicatesDimensions() {
    // Then, update the parameters of the ReplicatesCrystal and possibly the number of replicates.
    int ll = 1;
    int mm = 1;
    int nn = 1;

    while (unitCell.interfacialRadiusA * ll < cutOff) {
      ll++;
    }
    while (unitCell.interfacialRadiusB * mm < cutOff) {
      mm++;
    }
    while (unitCell.interfacialRadiusC * nn < cutOff) {
      nn++;
    }
    if (super.changeUnitCellParameters(
        unitCell.a * ll,
        unitCell.b * mm,
        unitCell.c * nn,
        unitCell.alpha,
        unitCell.beta,
        unitCell.gamma)) {
      l = ll;
      m = mm;
      n = nn;
      updateReplicateOperators();
      return true;
    }

    return false;
  }

  /**
   * Update the list of symmetry operators used to generate the replicates super cell from the
   * asymmetric unit.
   */
  private void updateReplicateOperators() {
    List<SymOp> symOps = spaceGroup.symOps;

    // First, we remove the existing symmetry operators.
    symOps.clear();

    /*
     Now create symmetry operators for each replicate. Note that the first
     symOp is still equivalent to the asymmetric unit and the first set of
     symOps are still equivalent to the unit cell.
    */
    double dX = 1.0 / (double) l;
    double dY = 1.0 / (double) m;
    double dZ = 1.0 / (double) n;
    for (int i = 0; i < l; i++) {
      for (int j = 0; j < m; j++) {
        for (int k = 0; k < n; k++) {
          int ii = 0;
          for (SymOp symOp : unitCell.spaceGroup.symOps) {
            double[] repTrans = new double[3];
            repTrans[0] = (symOp.tr[0] + i) * dX;
            repTrans[1] = (symOp.tr[1] + j) * dY;
            repTrans[2] = (symOp.tr[2] + k) * dZ;
            SymOp repSymOp = new SymOp(symOp.rot, repTrans);
            symOps.add(repSymOp);
            if (logger.isLoggable(Level.FINEST)) {
              logger.finest(format("\n SymOp (%2d,%2d,%2d): %d", i, j, k, ii));
              logger.finest(repSymOp.toString());
            }
            ii++;
          }
        }
      }
    }
  }
}
