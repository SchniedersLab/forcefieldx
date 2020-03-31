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
     * <p>
     * If there is no prior step (e.g. a newly loaded system that has not yet been rigidified), xPrior
     * and xNew can be copies of each other when passed to the method.
     * <p>
     * xPrior corresponds to atomCoordinates in OpenMM's constraint code. Ours will be in Angstroms, not nm.
     * xNew corresponds to atomCoordinatesP in OpenMM's constraint code. Ours will be in Angstroms, not nm.
     *
     * @param xPrior Atomic coordinates prior to the timestep to be constrained.
     * @param xNew   Atomic coordinates after the timestep; updated in-place to satisfy the constraint.
     * @param masses Masses.
     * @param tol    Acceptable constraint tolerance for numerical methods, as a fraction of bond length.
     */
    void applyConstraintToStep(final double[] xPrior, double[] xNew, final double[] masses, double tol);

    /**
     * Applies this Constraint to velocities, ensuring relative velocities are perpendicular to constrained
     * bonds, etc, without affecting positions. All arrays are globally indexed (i.e. includes all system
     * atoms, not just the constrained ones).
     * <p>
     * Our positions will be in Angstroms, and velocities in Angstroms/ps, compared to OpenMM's nm and nm/ps
     *
     * @param x      Atomic coordinates (unchanged).
     * @param v      Velocities (updated in-place to satisfy constraints).
     * @param masses Masses.
     * @param tol    Acceptable constraint tolerance for numerical methods; likely in Angstroms/ps
     */
    void applyConstraintToVelocities(final double[] x, double[] v, final double[] masses, double tol);

    /**
     * Returns the atomic XYZ indices of all Atoms constrained. Guaranteed to be unique.
     * The primary assumption will be that variables are in sets of 3x Cartesian coordinates.
     *
     * @return All indices of constrained Atoms.
     */
    int[] constrainedAtomIndices();

    /**
     * Checks if this Constraint is satisfied.
     *
     * @param x   Input coordinates to check.
     * @param tol Numerical tolerance as a fraction of bond stretch.
     * @return Whether this Constraint is satisfied.
     */
    boolean constraintSatisfied(final double[] x, double tol);

    /**
     * Checks if this Constraint is satisfied. Also checks velocities; bond
     * constraints, for example, require that relative velocity be orthogonal
     * to the bond. If the velocities vector is null or the tolerance is zero,
     * velocity checks are skipped.
     *
     * @param x    Input coordinates to check.
     * @param v    Input velocities to check. If null, velocity check disabled.
     * @param xTol Numerical tolerance for bond lengths.
     * @param vTol Numerical tolerance for velocity checks (typically in degrees). If zero, velocity check disabled.
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
