/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.utilities;

import java.lang.reflect.Array;

/**
 * Collected methods which allow easy implementation of <code>hashCode</code>.
 *
 * Example use case:
 * <pre>
 *  public int hashCode(){
 *    int result = HashCodeUtil.SEED;
 *    //collect the contributions of various fields
 *    result = HashCodeUtil.hash(result, fPrimitive);
 *    result = HashCodeUtil.hash(result, fObject);
 *    result = HashCodeUtil.hash(result, fArray);
 *    return result;
 *  }
 * </pre>
 *
 * @author Timothy D. Fenn
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class HashCodeUtil {

    /**
     * An initial value for a <code>hashCode</code>, to which is added
     * contributions from fields. Using a non-zero value decreases collisions of
     * <code>hashCode</code> values.
     */
    public static final int SEED = 23;

    /**
     * booleans.
     *
     * @param aSeed a int.
     * @param aBoolean a boolean.
     * @return a int.
     */
    public static int hash(int aSeed, boolean aBoolean) {
        return firstTerm(aSeed) + (aBoolean ? 1 : 0);
    }

    /**
     * chars.
     *
     * @param aSeed a int.
     * @param aChar a char.
     * @return a int.
     */
    public static int hash(int aSeed, char aChar) {
        return firstTerm(aSeed) + (int) aChar;
    }

    /**
     * ints.
     *
     * @param aSeed a int.
     * @param aInt a int.
     * @return a int.
     */
    public static int hash(int aSeed, int aInt) {
        /*
         * Implementation Note
         * Note that byte and short are handled by this method, through
         * implicit conversion.
         */
        return firstTerm(aSeed) + aInt;
    }

    /**
     * longs.
     *
     * @param aSeed a int.
     * @param aLong a long.
     * @return a int.
     */
    public static int hash(int aSeed, long aLong) {
        return firstTerm(aSeed) + (int) (aLong ^ (aLong >>> 32));
    }

    /**
     * floats.
     *
     * @param aSeed a int.
     * @param aFloat a float.
     * @return a int.
     */
    public static int hash(int aSeed, float aFloat) {
        return hash(aSeed, Float.floatToIntBits(aFloat));
    }

    /**
     * doubles.
     *
     * @param aSeed a int.
     * @param aDouble a double.
     * @return a int.
     */
    public static int hash(int aSeed, double aDouble) {
        return hash(aSeed, Double.doubleToLongBits(aDouble));
    }

    /**
     * <code>aObject</code> is a possibly-null object field, and possibly an
     * array.
     *
     * If <code>aObject</code> is an array, then each element may be a primitive
     * or a possibly-null object.
     *
     * @param aSeed a int.
     * @param aObject a {@link java.lang.Object} object.
     * @return a int.
     */
    public static int hash(int aSeed, Object aObject) {
        int result = aSeed;
        if (aObject == null) {
            result = hash(result, 0);
        } else if (!isArray(aObject)) {
            result = hash(result, aObject.hashCode());
        } else {
            int length = Array.getLength(aObject);
            for (int idx = 0; idx < length; ++idx) {
                Object item = Array.get(aObject, idx);
                //recursive call!
                result = hash(result, item);
            }
        }
        return result;
    }

    /// PRIVATE ///
    private static final int fODD_PRIME_NUMBER = 37;

    private static int firstTerm(int aSeed) {
        return fODD_PRIME_NUMBER * aSeed;
    }

    private static boolean isArray(Object aObject) {
        return aObject.getClass().isArray();
    }
}
