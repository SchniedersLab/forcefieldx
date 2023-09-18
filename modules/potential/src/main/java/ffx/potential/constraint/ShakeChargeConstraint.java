package ffx.potential.constraint;

import ffx.numerics.Constraint;

public class ShakeChargeConstraint implements Constraint {
    final int nConstraints;
    final double tol;
    final int c;

    public ShakeChargeConstraint(int nConstraints, int c, double tol){
        this.nConstraints = nConstraints;
        this.c = c;
        this.tol = tol;
    }
    @Override
    public void applyConstraintToStep(final double[] xPrior, double[] xNew, final double[] masses, double tol){
    }

    public boolean applyChargeConstraintToStep( final double[] xNew, final double[] accel, final double[] masses, final double dt){
        boolean done = true;

        double totalLambda = 0.0;
        double totalInverseMass = 0.0;
        for (int i = 0; i < nConstraints; i++) {
            double lambda = Math.sin(xNew[i]) * Math.sin(xNew[i]);
            totalLambda += lambda;
            totalInverseMass += (1.0 / masses[i]);
        }
        double delta = totalLambda - c;

        if (Math.abs(delta) > tol) {
            done = false;
            double term = delta / (dt * dt * totalInverseMass);
            for (int i = 0; i < nConstraints; i++) {
                accel[i] += -term * Math.sin(2 * xNew[i]);
            }
        }
        return done;
    }

    @Override
    public void applyConstraintToVelocities(final double[] x, double[] v, final double[] masses, double tol){

    }

    @Override
    public int[] constrainedAtomIndices() {
        return new int[0];
    }

    @Override
    public boolean constraintSatisfied(double[] x, double tol) {
        return false;
    }

    @Override
    public boolean constraintSatisfied(double[] x, double[] v, double xTol, double vTol) {
        return false;
    }

    @Override
    public int getNumDegreesFrozen() {
        return 0;
    }
}


