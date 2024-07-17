package ffx.numerics.multipole;

import java.util.HashMap;

import static ffx.numerics.math.ScalarMath.binomial;
import static java.lang.String.format;

public class MultipoleUtilities {


  /**
   * Returns the number of tensors for derivatives to the given order.
   *
   * @param order maximum number of derivatives.
   * @return the number of tensors.
   * @since 1.0
   */
  public static int tensorCount(int order) {
    long ret = binomial(order + 3, 3);
    assert (ret < Integer.MAX_VALUE);
    return (int) ret;
  }

  /**
   * Convenience method for writing out tensor indices.
   *
   * @param l number of d/dx partial derivatives.
   * @param m number of d/dx partial derivatives.
   * @param n number of d/dx partial derivatives.
   * @return a String of the form <code>Rlmn</code>.
   */
  protected static String rlmn(int l, int m, int n) {
    return format("R%d%d%d", l, m, n);
  }

  /**
   * Convenience method for writing out intermediate terms in the recursion.
   *
   * @param l number of d/dx partial derivatives.
   * @param m number of d/dx partial derivatives.
   * @param n number of d/dx partial derivatives.
   * @return a String of the form <code>termLMN</code>.
   */
  protected static String term(int l, int m, int n) {
    return format("term%d%d%d", l, m, n);
  }

  /**
   * Convenience method for writing out intermediate terms in the recursion.
   *
   * @param l number of d/dx partial derivatives.
   * @param m number of d/dx partial derivatives.
   * @param n number of d/dx partial derivatives.
   * @param j the jth intermediate term.
   * @return a String of the form <code>termLMNJ</code>.
   */
  protected static String term(int l, int m, int n, int j) {
    return format("term%d%d%d%d", l, m, n, j);
  }

  /**
   * The index is based on the idea of filling tetrahedron.
   * <p>
   * 1/r has an index of 0.
   * <br>
   * derivatives of x are first; indices from 1..o for d/dx..(d/dx)^o
   * <br>
   * derivatives of x and y are second; base triangle of size (o+1)(o+2)/2
   * <br>
   * derivatives of x, y and z are last; total size (o+1)*(o+2)*(o+3)/6
   * <br>
   * <p>
   * This function is useful to set up masking constants:
   * <br>
   * static int Tlmn = ti(l,m,n,order)
   * <br>
   * For example the (d/dy)^2 (1/R) storage location:
   * <br>
   * static int T020 = ti(0,2,0,order)
   *
   * @param dx    int The number of d/dx operations.
   * @param dy    int The number of d/dy operations.
   * @param dz    int The number of d/dz operations.
   * @param order int The maximum tensor order (0 .LE. dx + dy + dz .LE. order).
   * @return int in the range (0..binomial(order + 3, 3) - 1)
   */
  protected static int ti(int dx, int dy, int dz, int order) {
    if (dx < 0 || dy < 0 || dz < 0 || dx + dy + dz > order) {
      return -1;
    }

    int size = (order + 1) * (order + 2) * (order + 3) / 6;
    /*
     We only get to the top of the tetrahedron if dz = order, otherwise
     subtract off the top, including the level of the requested tensor
     index.
    */
    int top = order + 1 - dz;
    top = top * (top + 1) * (top + 2) / 6;
    int zIndex = size - top;
    /*
     Given the "dz level", dy can range from (0 .. order - dz).
     To get to the row for a specific value of dy, dy*(order + 1) - dy*(dy-1)/2 indices are skipped.
     This is an operation that looks like the area of rectangle, minus the area of an empty triangle.
    */
    int yIndex = dy * (order - dz) - (dy - 1) * (dy - 2) / 2 + 1;
    /*
     Given the dz level and dy row, dx can range from (0..order - dz - dy)
     The dx index is just walking down the dy row for "dx" steps.
    */
    return dx + yIndex + zIndex;
  }

  /**
   * Convenience method for writing out tensor indices.
   *
   * @param l number of d/dx partial derivatives.
   * @param m number of d/dx partial derivatives.
   * @param n number of d/dx partial derivatives.
   * @return a String of the form <code>lmn</code>.
   */
  protected static String lmn(int l, int m, int n) {
    return format("%d%d%d", l, m, n);
  }

  /**
   * Load a tensor element into a SIMD register.
   * <p>
   * For l=m=n=0, the returned String is:
   * DoubleVector t000 = fromArray(SPECIES, t, T000);
   *
   * @param l The number of d/dx operations.
   * @param n The number of d/dy operations.
   * @param m The number of d/dz operations.
   * @return a String of the form <code>DoubleVector t000 = fromArray(SPECIES, t, T000);</code>
   */
  protected static String loadTensor(int l, int n, int m, HashMap<Integer, String> tensorMap) {
    String s = lmn(l, n, m);
    String offset = "T" + lmn(l, n, m);
    String name = "t" + s;
    if (tensorMap.containsValue(name)) {
      return "";
    }
    int size = tensorMap.size();
    tensorMap.put(size, name);
    return format("\tDoubleVector %s = fromArray(SPECIES, t, %s);\n", name, offset);
  }

  /**
   * Code to store an electrostatic potential element into an array.
   * <p>
   * For l=m=n=0, the returned String is (including a preceding tab):
   * term000.intoArray(to, T000);
   *
   * @param to variable name of the array to store into.
   * @param l  number of d/dx partial derivatives.
   * @param n  number of d/dy partial derivatives.
   * @param m  number of d/dz partial derivatives.
   * @return a String of the form <code>term000.intoArray(to, T000);</code>
   */
  protected static String storePotential(String to, int l, int n, int m) {
    String from = term(l, n, m);
    String offset = "T" + lmn(l, n, m);
    return format("\t%s.intoArray(%s, %s);\n", from, to, offset);
  }

  /**
   * Code to store a negated electrostatic potential element into an array.
   * <p>
   * For l=m=n=0, the returned String is (including a preceding tab):
   * term000.neg().intoArray(to, T000);
   *
   * @param to variable name of the array to store into.
   * @param l  number of d/dx partial derivatives.
   * @param n  number of d/dy partial derivatives.
   * @param m  number of d/dz partial derivatives.
   * @return a String of the form <code>term000.intoArray(to, T000);</code>
   */
  protected static String storePotentialNeg(String to, int l, int n, int m) {
    String from = term(l, n, m);
    String offset = "T" + lmn(l, n, m);
    return format("\t%s.neg().intoArray(%s, %s);\n", from, to, offset);
  }
}
