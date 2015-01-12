/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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
package ffx.numerics;

/**
 * <p>
 * Scalar class.</p>
 *
 * @author Timothy D. Fenn
 *
 * @since 1.0
 */
public class Scalar {

    private final double xyz[];
    private double scalar;

    /**
     * <p>
     * Constructor for Scalar.</p>
     */
    public Scalar() {
        xyz = new double[3];
    }

    /**
     * <p>
     * Constructor for Scalar.</p>
     *
     * @param xyz an array of double.
     * @param scalar a double.
     */
    public Scalar(double xyz[], double scalar) {
        this(xyz[0], xyz[1], xyz[2], scalar);
    }

    /**
     * <p>
     * Constructor for Scalar.</p>
     *
     * @param xyz an array of int.
     * @param scalar a double.
     */
    public Scalar(int xyz[], double scalar) {
        this((double) xyz[0], (double) xyz[1], (double) xyz[2], scalar);
    }

    /**
     * <p>
     * Constructor for Scalar.</p>
     *
     * @param x a int.
     * @param y a int.
     * @param z a int.
     * @param scalar a double.
     */
    public Scalar(int x, int y, int z, double scalar) {
        this((double) x, (double) y, (double) z, scalar);
    }

    /**
     * <p>
     * Constructor for Scalar.</p>
     *
     * @param x a double.
     * @param y a double.
     * @param z a double.
     * @param scalar a double.
     */
    public Scalar(double x, double y, double z, double scalar) {
        this();
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
        this.scalar = scalar;
    }

    /**
     * <p>
     * getXYZ</p>
     *
     * @return an array of double.
     */
    public double[] getXYZ() {
        return xyz;
    }

    /**
     * <p>
     * setXYZ</p>
     *
     * @param xyz an array of double.
     */
    public void setXYZ(double xyz[]) {
        this.xyz[0] = xyz[0];
        this.xyz[1] = xyz[1];
        this.xyz[2] = xyz[2];
    }

    /**
     * <p>
     * Getter for the field <code>scalar</code>.</p>
     *
     * @return a double.
     */
    public double getScalar() {
        return scalar;
    }

    /**
     * <p>
     * Setter for the field <code>scalar</code>.</p>
     *
     * @param scalar a double.
     */
    public void setScalar(double scalar) {
        this.scalar = scalar;
    }
}
