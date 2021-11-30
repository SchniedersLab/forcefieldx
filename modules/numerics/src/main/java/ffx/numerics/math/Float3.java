// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
// ******************************************************************************
package ffx.numerics.math;

import static java.lang.System.arraycopy;

/**
 * Convenience class for working with 3D float vectors.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Float3 {

  /** Internal storage for the vector. */
  private final float[] a;

  /** Construct a Float3 at (0.0, 0.0, 0.0). */
  public Float3() {
    a = new float[] {0.0f, 0.0f, 0.0f};
  }

  /**
   * Construct a Float3 at (x, y, z).
   *
   * @param x X value.
   * @param y Y value.
   * @param z Z value.
   */
  public Float3(float x, float y, float z) {
    a = new float[] {x, y, z};
  }

  /**
   * Construct a Float3 at a.
   *
   * @param a Float3 to initialize from (the data is copied).
   */
  public Float3(float[] a) {
    this.a = new float[] {a[0], a[1], a[2]};
  }

  /**
   * Cross product of this Float3 with b.
   *
   * @param b Vector b.
   * @return Returns the cross product in a new Float3.
   */
  public Float3 X(Float3 b) {
    return new Float3(FloatMath.X(a, b.get()));
  }

  /**
   * In-place Cross product of this Float3 with b.
   *
   * @param b Vector b.
   * @return Returns the cross product in this Float3.
   */
  public Float3 XI(Float3 b) {
    return set(X(b));
  }

  /**
   * Finds the sum of this Float3 with b.
   *
   * @param b Second Float3.
   * @return Returns the sum in a new Float3.
   */
  public Float3 add(Float3 b) {
    return new Float3(FloatMath.add(a, b.get()));
  }

  /**
   * Finds the sum of this Float3 with b in place.
   *
   * @param b Second Float3.
   * @return Returns the sum in this Float3.
   */
  public Float3 addI(Float3 b) {
    return set(add(b));
  }

  /**
   * Angle of this Float3 with b.
   *
   * @param b Vector b.
   * @return Returns the angle.
   */
  public float angle(Float3 b) {
    return FloatMath.angle(a, b.get());
  }

  /**
   * Returns a new copy of this Float3.
   *
   * @return Returns a copy of this Float3.
   */
  public Float3 copy() {
    return new Float3(a);
  }

  /**
   * Finds the distance between two vectors.
   *
   * @param b Second vector.
   * @return Returns the distance between this Float3 and b.
   */
  public float dist(Float3 b) {
    return FloatMath.dist(a, b.get());
  }

  /**
   * Finds the squared distance between two vectors
   *
   * @param b Second vector.
   * @return Returns the squared distance between this Float3 and b.
   */
  public float dist2(Float3 b) {
    return FloatMath.dist2(a, b.get());
  }

  /**
   * Finds the dot product between two vectors.
   *
   * @param b Second vector.
   * @return Returns the dot product of this Float3 and b.
   */
  public float dot(Float3 b) {
    return FloatMath.dot(a, b.get());
  }

  /**
   * Returns a reference to the internal float array that stores this Float3.
   *
   * @return A reference to the internal float array.
   */
  public float[] get() {
    return a;
  }

  /**
   * Finds the length of this Float3.
   *
   * @return Length of vector this Float3.
   */
  public float length() {
    return FloatMath.length(a);
  }

  /**
   * Finds the length of this Float3 squared.
   *
   * @return Length of vector this Float3 squared.
   */
  public float length2() {
    return FloatMath.length2(a);
  }

  /** Log this Float3. */
  public void log() {
    FloatMath.log(a);
  }

  /**
   * Normalize this Float3.
   *
   * @return Returns a normalized Float3.
   */
  public Float3 normalize() {
    return new Float3(FloatMath.normalize(a));
  }

  /**
   * Normalize this Float3 in place.
   *
   * @return Returns a reference to this Float3 normalized.
   */
  public Float3 normalizeI() {
    return set(normalize());
  }

  /**
   * Scales a Float3.
   *
   * @param d A scalar value.
   * @return Returns a new scaled Float3.
   */
  public Float3 scale(float d) {
    return new Float3(FloatMath.scale(a, d));
  }

  /**
   * Scales a Float3 in place.
   *
   * @param d A scalar value.
   * @return Returns a reference to this Float3 scaled.
   */
  public Float3 scaleI(float d) {
    return set(scale(d));
  }

  /**
   * Set the value of this Float3.
   *
   * @param x X-value.
   * @param y Y-value.
   * @param z Z-value.
   * @return A reference to this Float3.
   */
  public Float3 set(float x, float y, float z) {
    a[0] = x;
    a[1] = y;
    a[2] = z;
    return this;
  }

  /**
   * Set the value of this Float3.
   *
   * @param b Double array that is copied.
   * @return A reference to this Float3.
   */
  public Float3 set(float[] b) {
    arraycopy(b, 0, a, 0, 3);
    return this;
  }

  /**
   * Set the value of this Float3.
   *
   * @param b Float3 that is copied.
   * @return A reference to this Float3.
   */
  public Float3 set(Float3 b) {
    arraycopy(b.get(), 0, a, 0, 3);
    return this;
  }

  /**
   * Finds the difference between two vectors.
   *
   * @param b Second vector
   * @return Returns the difference in a new Float3.
   */
  public Float3 sub(Float3 b) {
    return new Float3(FloatMath.sub(a, b.get()));
  }

  /**
   * Finds the difference between two vectors.
   *
   * @param b Second vector
   * @return Returns the difference in this Float3.
   */
  public Float3 subI(Float3 b) {
    return set(sub(b));
  }

  /**
   * Describe this Float3 in a String.
   *
   * @return Returns a String description.
   */
  @Override
  public String toString() {
    return FloatMath.toString(a);
  }

  /**
   * Get the value x.
   *
   * @return Returns x.
   */
  public float x() {
    return a[0];
  }

  /**
   * Get the value of y.
   *
   * @return Returns y.
   */
  public float y() {
    return a[1];
  }

  /**
   * Get the value of z.
   *
   * @return Returns z.
   */
  public float z() {
    return a[2];
  }
}
