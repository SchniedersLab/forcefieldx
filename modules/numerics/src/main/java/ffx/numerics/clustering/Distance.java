// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.numerics.clustering;

/**
 * Simple value object storing a distance and an optional weight used during
 * linkage computations; comparable by distance and cloneable.
 *
 * @author Lars Behnke, 2013
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Distance implements Comparable<Distance>, Cloneable {

  private Double distance;
  private Double weight;

  /**
   * Creates a Distance with value 0 and weight 1.
   */
  public Distance() {
    this(0.0);
  }

  /**
   * Creates a Distance with the given value and unit weight.
   *
   * @param distance the distance value
   */
  public Distance(Double distance) {
    this(distance, 1.0);
  }

  /**
   * Creates a Distance with the given value and weight.
   *
   * @param distance the distance value
   * @param weight   the weight associated with this distance
   */
  public Distance(Double distance, Double weight) {
    this.distance = distance;
    this.weight = weight;
  }

  /**
   * Gets the distance value.
   *
   * @return the distance value (may be null)
   */
  public Double getDistance() {
    return distance;
  }

  /**
   * Sets the distance value.
   *
   * @param distance the distance value to set
   */
  public void setDistance(Double distance) {
    this.distance = distance;
  }

  /**
   * Gets the weight.
   *
   * @return the weight
   */
  public Double getWeight() {
    return weight;
  }

  /**
   * Sets the weight.
   *
   * @param weight the weight to set
   */
  public void setWeight(Double weight) {
    this.weight = weight;
  }

  /**
   * Checks whether the distance value is NaN or undefined.
   *
   * @return true if the distance is null or NaN; false otherwise
   */
  public boolean isNaN() {
    return distance == null || distance.isNaN();
  }

  /**
   * Compares by distance value, with null other treated as greater (this &lt; other).
   *
   * @param distance the other Distance to compare to
   * @return negative, zero, or positive per Comparable contract
   */
  @Override
  public int compareTo(Distance distance) {
    return distance == null ? 1 : getDistance().compareTo(distance.getDistance());
  }

  /**
   * Returns a string containing distance and weight values.
   *
   * @return formatted string representation
   */
  @Override
  public String toString() {
    return String.format("distance : %.2f, weight : %.2f", distance, weight);
  }

  /**
   * Creates a clone of this Distance.
   *
   * @return a clone of this Distance
   */
  @Override
  public Distance clone() {
    try {
      Distance clone = (Distance) super.clone();
      clone.distance = distance;
      clone.weight = weight;
      return clone;
    } catch (CloneNotSupportedException e) {
      throw new AssertionError();
    }
  }
}
