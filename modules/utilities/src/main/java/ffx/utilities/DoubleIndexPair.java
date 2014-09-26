/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.utilities;

/**
 * <p>
 * DoubleIndexPair class.</p>
 *
 * @author Jacob M. Litman
 *
 * @since 1.0
 */
public class DoubleIndexPair implements Comparable {

    private final int index;
    private final double doubleValue;

    /**
     * Allows sorting of floating-point values while retaining knowledge of
     * where that value was in some original list or array.
     *
     * @param index
     * @param doubleValue
     */
    public DoubleIndexPair(int index, double doubleValue) {
        this.index = index;
        this.doubleValue = doubleValue;
    }

    public int getIndex() {
        return index;
    }

    public double getDoubleValue() {
        return doubleValue;
    }

    @Override
    public int compareTo(Object o) {
        if (o == null) {
            return 0;
        }
        if (!(o instanceof DoubleIndexPair)) {
            return 0;
        }
        DoubleIndexPair other = (DoubleIndexPair) o;
        if (doubleValue < other.doubleValue) {
            return -1;
        } else if (doubleValue > other.doubleValue) {
            return 1;
        } else {
            return 0;
        }
    }
}
