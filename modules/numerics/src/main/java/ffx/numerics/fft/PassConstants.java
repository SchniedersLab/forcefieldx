package ffx.numerics.fft;

/**
 * Constant factors needed for each pass.
 *
 * @param factor         The factor.
 * @param product        The product of all factors applied so far.
 * @param outerLoopLimit The outer loop limit (n / product).
 * @param innerLoopLimit The inner loop limit (product / factor).
 * @param nextInput      The next input (n / factor).
 * @param di             Twice the next input to account for complex numbers.
 * @param dj             Twice the inner loop limit to account for complex numbers.
 * @param twiddles       The twiddle factors for this pass.
 */
public record PassConstants(int factor, int product, int outerLoopLimit, int innerLoopLimit, int nextInput,
                            int di, int dj, double[][] twiddles) {
  // Empty.
}
