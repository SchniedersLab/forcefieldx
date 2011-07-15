/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT Aexty WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.numerics;

/**
 *
 * @author Tim Fenn
 */
public class Scalar {

    private double xyz[] = new double[3];
    private double scalar;

    public Scalar() {
    }

    public Scalar(double xyz[], double scalar) {
        this(xyz[0], xyz[1], xyz[2], scalar);
    }

    public Scalar(int xyz[], double scalar) {
        this((double) xyz[0], (double) xyz[1], (double) xyz[2], scalar);
    }

    public Scalar(int x, int y, int z, double scalar) {
        this((double) x, (double) y, (double) z, scalar);
    }

    public Scalar(double x, double y, double z, double scalar) {
        this.xyz[0] = x;
        this.xyz[1] = y;
        this.xyz[2] = z;
        this.scalar = scalar;
    }

    public double[] getXYZ() {
        return xyz;
    }

    public void setXYZ(double xyz[]) {
        this.xyz[0] = xyz[0];
        this.xyz[1] = xyz[1];
        this.xyz[2] = xyz[2];
    }

    public double getScalar() {
        return scalar;
    }

    public void setScalar(double scalar) {
        this.scalar = scalar;
    }
}
