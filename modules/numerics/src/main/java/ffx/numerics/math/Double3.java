//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.numerics.math;

import static java.lang.System.arraycopy;

/**
 * Convenience class for working with 3D double vectors.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Double3 {

    /**
     * Internal storage for the vector.
     */
    private final double[] a;

    /**
     * Construct a Double3 at (0.0, 0.0, 0.0).
     */
    public Double3() {
        a = new double[]{0.0, 0.0, 0.0};
    }

    /**
     * Construct a Double3 at (x, y, z).
     *
     * @param x X value.
     * @param y Y value.
     * @param z Z value.
     */
    public Double3(double x, double y, double z) {
        a = new double[]{x, y, z};
    }

    /**
     * Construct a Double3 at a.
     *
     * @param a Double3 to initialize from (the data is copied).
     */
    public Double3(double[] a) {
        this.a = new double[]{a[0], a[1], a[2]};
    }

    /**
     * Cross product of this Double3 with b.
     *
     * @param b Vector b.
     * @return Returns the cross product in a new Double3.
     */
    public Double3 X(Double3 b) {
        return new Double3(DoubleMath.X(a, b.get()));
    }

    /**
     * In-place Cross product of this Double3 with b.
     *
     * @param b Vector b.
     * @return Returns the cross product in this Double3.
     */
    public Double3 XI(Double3 b) {
        return set(X(b));
    }

    /**
     * Finds the sum of this Double3 with b.
     *
     * @param b Second Double3.
     * @return Returns the sum in a new Double3.
     */
    public Double3 add(Double3 b) {
        return new Double3(DoubleMath.add(a, b.get()));
    }

    /**
     * Finds the sum of this Double3 with b in place.
     *
     * @param b Second Double3.
     * @return Returns the sum in this Double3.
     */
    public Double3 addI(Double3 b) {
        return set(add(b));
    }

    /**
     * Angle of this Double3 with b.
     *
     * @param b Vector b.
     * @return Returns the angle.
     */
    public double angle(Double3 b) {
        return DoubleMath.angle(a, b.get());
    }

    /**
     * Returns a new copy of this Double3.
     *
     * @return Returns a copy of this Double3.
     */
    public Double3 copy() {
        return new Double3(a);
    }

    /**
     * Finds the Euclidean distance between two positions.
     *
     * @param b Second vector.
     * @return Returns the distance between this Double3 and b.
     */
    public double dist(Double3 b) {
        return DoubleMath.dist(a, b.get());
    }

    /**
     * Finds the square of the Euclidean distance between two postions.
     *
     * @param b Second vector.
     * @return Returns the squared distance between this Double3 and b.
     */
    public double dist2(Double3 b) {
        return DoubleMath.dist2(a, b.get());
    }

    /**
     * Finds the dot product between two vectors.
     *
     * @param b Second vector.
     * @return Returns the dot product of this Double3 and b.
     */
    public double dot(Double3 b) {
        return DoubleMath.dot(a, b.get());
    }

    /**
     * Returns a reference to the internal double array that stores this Double3.
     *
     * @return A reference to the internal double array.
     */
    public double[] get() {
        return a;
    }

    /**
     * Finds the length of this Double3.
     *
     * @return Length of vector this Double3.
     */
    public double length() {
        return DoubleMath.length(a);
    }

    /**
     * Finds the length of this Double3 squared.
     *
     * @return Length of vector this Double3 squared.
     */
    public double length2() {
        return DoubleMath.length2(a);
    }

    /**
     * Log this Double3.
     */
    public void log() {
        DoubleMath.log(a);
    }

    /**
     * Normalize this Double3.
     *
     * @return Returns a normalized Double3.
     */
    public Double3 normalize() {
        return new Double3(DoubleMath.normalize(a));
    }

    /**
     * Normalize this Double3 in place.
     *
     * @return Returns a reference to this Double3 normalized.
     */
    public Double3 normalizeI() {
        return set(normalize());
    }

    /**
     * Scales a Double3.
     *
     * @param d A scalar value.
     * @return Returns a new scaled Double3.
     */
    public Double3 scale(double d) {
        return new Double3(DoubleMath.scale(a, d));
    }

    /**
     * Scales a Double3 in place.
     *
     * @param d A scalar value.
     * @return Returns a reference to this Double3 scaled.
     */
    public Double3 scaleI(double d) {
        return set(scale(d));
    }

    /**
     * Set the value of this Double3.
     *
     * @param x X-value.
     * @param y Y-value.
     * @param z Z-value.
     * @return A reference to this Double3.
     */
    public Double3 set(double x, double y, double z) {
        a[0] = x;
        a[1] = y;
        a[2] = z;
        return this;
    }

    /**
     * Set the value of this Double3.
     *
     * @param b Double array that is copied.
     * @return A reference to this Double3.
     */
    public Double3 set(double[] b) {
        arraycopy(b, 0, a, 0, 3);
        return this;
    }

    /**
     * Set the value of this Double3.
     *
     * @param b Double3 that is copied.
     * @return A reference to this Double3.
     */
    public Double3 set(Double3 b) {
        arraycopy(b.get(), 0, a, 0, 3);
        return this;
    }

    /**
     * Finds the difference between two vectors.
     *
     * @param b Second vector
     * @return Returns the difference in a new Double3.
     */
    public Double3 sub(Double3 b) {
        return new Double3(DoubleMath.sub(a, b.get()));
    }

    /**
     * Finds the difference between two vectors.
     *
     * @param b Second vector
     * @return Returns the difference in this Double3.
     */
    public Double3 subI(Double3 b) {
        return set(sub(b));
    }

    /**
     * Describe this Double3 in a String.
     *
     * @return Returns a String description.
     */
    @Override
    public String toString() {
        return DoubleMath.toString(a);
    }

    /**
     * Get the value x.
     *
     * @return Returns x.
     */
    public double x() {
        return a[0];
    }

    /**
     * Get the value of y.
     *
     * @return Returns y.
     */
    public double y() {
        return a[1];
    }

    /**
     * Get the value of z.
     *
     * @return Returns z.
     */
    public double z() {
        return a[2];
    }
}
    
    
