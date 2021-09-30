//******************************************************************************
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
//******************************************************************************
package ffx.numerics;

/**
 * Defines a set of geometric constraints that must be applied self-consistently.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface Constraint {

  /**
   * Applies this Constraint in the context of a partially calculated MD timestep. All arrays are
   * globally indexed (i.e. includes all system atoms, not just the constrained ones).
   *
   * <p>If there is no prior step (e.g. a newly loaded system that has not yet been rigidified),
   * xPrior and xNew can be copies of each other when passed to the method.
   *
   * <p>xPrior corresponds to atomCoordinates in OpenMM's constraint code. Ours will be in
   * Angstroms, not nm. xNew corresponds to atomCoordinatesP in OpenMM's constraint code. Ours will
   * be in Angstroms, not nm.
   *
   * @param xPrior Atomic coordinates prior to the timestep to be constrained.
   * @param xNew Atomic coordinates after the timestep; updated in-place to satisfy the constraint.
   * @param masses Masses.
   * @param tol Acceptable constraint tolerance for numerical methods, as a fraction of bond length.
   */
  void applyConstraintToStep(
      final double[] xPrior, double[] xNew, final double[] masses, double tol);

  /**
   * Applies this Constraint to velocities, ensuring relative velocities are perpendicular to
   * constrained bonds, etc, without affecting positions. All arrays are globally indexed (i.e.
   * includes all system atoms, not just the constrained ones).
   *
   * <p>Our positions will be in Angstroms, and velocities in Angstroms/ps, compared to OpenMM's nm
   * and nm/ps
   *
   * @param x Atomic coordinates (unchanged).
   * @param v Velocities (updated in-place to satisfy constraints).
   * @param masses Masses.
   * @param tol Acceptable constraint tolerance for numerical methods; likely in Angstroms/ps
   */
  void applyConstraintToVelocities(final double[] x, double[] v, final double[] masses, double tol);

  /**
   * Returns the atomic XYZ indices of all Atoms constrained. Guaranteed to be unique. The primary
   * assumption will be that variables are in sets of 3x Cartesian coordinates.
   *
   * @return All indices of constrained Atoms.
   */
  int[] constrainedAtomIndices();

  /**
   * Checks if this Constraint is satisfied.
   *
   * @param x Input coordinates to check.
   * @param tol Numerical tolerance as a fraction of bond stretch.
   * @return Whether this Constraint is satisfied.
   */
  boolean constraintSatisfied(final double[] x, double tol);

  /**
   * Checks if this Constraint is satisfied. Also checks velocities; bond constraints, for example,
   * require that relative velocity be orthogonal to the bond. If the velocities vector is null or
   * the tolerance is zero, velocity checks are skipped.
   *
   * @param x Input coordinates to check.
   * @param v Input velocities to check. If null, velocity check disabled.
   * @param xTol Numerical tolerance for bond lengths.
   * @param vTol Numerical tolerance for velocity checks (typically in degrees). If zero, velocity
   *     check disabled.
   * @return Whether this Constraint is satisfied.
   */
  boolean constraintSatisfied(final double[] x, final double[] v, double xTol, double vTol);

  /**
   * Returns the number of degrees of freedom this Constraint constrains.
   *
   * @return Number of frozen DoF.
   */
  int getNumDegreesFrozen();
}
