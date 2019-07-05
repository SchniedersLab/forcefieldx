package ffx.potential.constraint;

import ffx.potential.bonded.Bond;

import java.util.List;

/**
 * Defines a set of geometric constraints that must be applied self-consistently.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface Constraint {
    /**
     * Returns the atomic XYZ indices of all Atoms constrained. Guaranteed to be unique.
     * @return All indices of constrained Atoms.
     */
    int[] constrainedAtomIndices();

    /**
     * Applies this Constraint in the context of a partially calculated MD timestep. All arrays are
     * globally indexed (i.e. includes all system atoms, not just the constrained ones).
     *
     * If there is no prior step (e.g. a newly loaded system that has not yet been rigidified), xPrior
     * and xNew can be copies of each other when passed to the method.
     *
     * xPrior corresponds to atomCoordinates in OpenMM's constraint code. Ours will be in Angstroms, not nm.
     * xNew corresponds to atomCoordinatesP in OpenMM's constraint code. Ours will be in Angstroms, not nm.
     *
     * @param xPrior Atomic coordinates prior to the timestep to be constrained.
     * @param xNew   Atomic coordinates after the timestep; updated in-place to satisfy the constraint.
     * @param masses Masses.
     * @param tol    Acceptable constraint tolerance for numerical methods, as a fraction of bond length.
     */
    void applyConstraintToStep(final double[] xPrior, double[] xNew, final double[] masses, double tol);

    // OpenMM also has an applyToVelocities method, useful for methods which calculate both a half-step and full-step
    // velocity (e.g. Velocity Verlet), so as to ensure velocities are constrained at both half-step and full-step.

    /**
     * Applies this Constraint to velocities, ensuring relative velocities are perpendicular to constrained
     * bonds, etc, without affecting positions. All arrays are globally indexed (i.e. includes all system
     * atoms, not just the constrained ones).
     *
     * Our positions will be in Angstroms, and velocities in Angstroms/ps, compared to OpenMM's nm and nm/ps
     *
     * @param x      Atomic coordinates (unchanged).
     * @param v      Velocities (updated in-place to satisfy constraints).
     * @param masses Masses.
     * @param tol    Acceptable constraint tolerance for numerical methods; likely in Angstroms/ps
     */
    void applyConstraintToVelocities(final double[] x, double[] v, final double[] masses, double tol);
}
