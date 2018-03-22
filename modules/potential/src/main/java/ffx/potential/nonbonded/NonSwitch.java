package ffx.potential.nonbonded;

/**
 * Implements a "dummy switch" for when no cutoff is desired in classes that depend on having a switching function to
 * smoothly go to 0. Not intended for DualTopologyEnergy: used primarily in VanDerWaals to replace MultiplicativeSwitch.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 */
public class NonSwitch extends MultiplicativeSwitch {

    /**
     * Constructs a non-switch where f(x) = 1 for all x.
     */
    public NonSwitch() {

    }

    /**
     * Value of a non-switch is always 1.
     * @param r Ignored.
     * @return Always 1
     */
    public double taper(double r) {
        return 1;
    }

    /**
     * Derivatives of a non-switch are always zero.
     * @param r Ignored.
     * @return Always 0
     */
    public double dtaper(double r) {
        return 0;
    }

    /**
     * Value of a non-switch is always 1.
     * @param r Ignored.
     * @param r2 Ignored.
     * @param r3 Ignored.
     * @param r4 Ignored.
     * @param r5 Ignored.
     * @return Always 1
     */
    public double taper(double r, double r2, double r3, double r4, double r5) {
        return 1;
    }

    /**
     * Derivatives of a non-switch are always zero.
     * @param r Ignored.
     * @param r2 Ignored.
     * @param r3 Ignored.
     * @param r4 Ignored.
     * @return Always 0.
     */
    public double dtaper(double r, double r2, double r3, double r4) {
        return 0;
    }

    @Override
    public double getZeroBound() {
        throw new UnsupportedOperationException(" This is not actually a switching function!");
    }

    @Override
    public double getOneBound() {
        throw new UnsupportedOperationException(" This is not actually a switching function!");
    }

    @Override
    public boolean constantOutsideBounds() {
        throw new UnsupportedOperationException(" This is not actually a switching function!");
    }

    @Override
    public boolean validOutsideBounds() {
        throw new UnsupportedOperationException(" This is not actually a switching function!");
    }

    @Override
    public int getHighestOrderZeroDerivative() {
        throw new UnsupportedOperationException(" This is not actually a switching function!");
    }

    @Override
    public boolean symmetricToUnity() {
        throw new UnsupportedOperationException(" This is not actually a switching function!");
    }

    /**
     * Value of a non-switch is always 1.
     *
     * @param x Ignored.
     * @return Always 0.
     * @throws IllegalArgumentException Never.
     */
    @Override
    public double valueAt(double x) throws IllegalArgumentException {
        return 1;
    }

    /**
     * Derivatives of a non-switch are always zero.
     * @param x Ignored.
     * @return Always 0.
     */
    @Override
    public double firstDerivative(double x) {
        return 0;
    }

    /**
     * Derivatives of a non-switch are always zero.
     * @param x Ignored.
     * @return Always 0.
     */
    @Override
    public double secondDerivative(double x) {
        return 0;
    }

    /**
     * Derivatives of a non-switch are always zero.
     * @param x Ignored.
     * @param order Ignored.
     * @return Always 0.
     */
    @Override
    public double nthDerivative(double x, int order) throws IllegalArgumentException {
        return 0;
    }

    @Override
    public String toString() {
        return "Constant, 1-valued non-switch that always returns 1 (for when no cutoff is desired).";
    }
}
