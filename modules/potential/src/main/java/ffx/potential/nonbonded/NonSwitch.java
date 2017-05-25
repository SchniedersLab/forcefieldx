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
     * Value of the switching function at r.
     * @param r r
     * @return Always 1
     */
    public double taper(double r) {
        return 1;
    }

    /**
     * First derivative of the switching function at r.
     * @param r r
     * @return Always 0
     */
    public double dtaper(double r) {
        return 0;
    }

    /**
     * Value of the switching function at r.
     * @param r r
     * @param r2 r^2
     * @param r3 r^3
     * @param r4 r^4
     * @param r5 r^5
     * @return Always 1
     */
    public double taper(double r, double r2, double r3, double r4, double r5) {
        return 1;
    }

    /**
     * First derivative of the switching function at r.
     * @param r r
     * @param r2 r^2
     * @param r3 r^3
     * @param r4 r^4
     * @return First derivative of switch at r
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

    @Override
    public double valueAt(double x) throws IllegalArgumentException {
        return 1;
    }

    @Override
    public double firstDerivative(double x) {
        return 0;
    }

    @Override
    public double secondDerivative(double x) {
        return 0;
    }

    @Override
    public double nthDerivative(double x, int order) throws IllegalArgumentException {
        return 0;
    }

    @Override
    public String toString() {
        return "Constant, 1-valued non-switch that always returns 1 (for when no cutoff is desired).";
    }
}
