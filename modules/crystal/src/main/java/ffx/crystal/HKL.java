// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
package ffx.crystal;

import java.util.Objects;

import static org.apache.commons.math3.util.FastMath.PI;

/**
 * The HKL class represents a single reflection.
 *
 * @author Timothy D. Fenn
 * @see ReflectionList
 * @since 1.0
 */
public class HKL {

  /**
   * Constant <code>ndiv=12.0</code>
   */
  static final double ndiv = 12.0;

  /**
   * The h-index of the reflection.
   */
  protected int h;
  /**
   * The k-index of the reflection.
   */
  protected int k;
  /**
   * The l-index of the reflection.
   */
  protected int l;
  /**
   * The epsilon value of the reflection, which is used for systematic absences.
   */
  protected int epsilon;
  /**
   * Allowed is used for centric reflections.
   */
  protected int allowed;
  /**
   * The bin number of this reflection, which is used for resolution dependent R/Rfree.
   */
  protected int bin;
  /**
   * The unique index of this reflection.
   */
  protected int index;

  /**
   * Constructor for HKL.
   */
  public HKL() {
  }

  /**
   * Constructor for HKL.
   *
   * @param h The h-index of the reflection.
   * @param k The k-index of the reflection.
   * @param l The l-index of the reflection.
   */
  public HKL(int h, int k, int l) {
    this.h = h;
    this.k = k;
    this.l = l;
    allowed = 255;
  }

  /**
   * Constructor for HKL.
   *
   * @param h       The h-index of the reflection.
   * @param k       The k-index of the reflection.
   * @param l       The l-index of the reflection.
   * @param eps     an int.
   * @param allowed an int.
   */
  public HKL(int h, int k, int l, int eps, int allowed) {
    this.h = h;
    this.k = k;
    this.l = l;
    this.epsilon = eps;
    this.allowed = allowed;
  }

  /**
   * Negate the reflection.
   *
   * @return a {@link ffx.crystal.HKL} object.
   */
  public HKL neg() {
    return new HKL(-h, -k, -l);
  }

  /**
   * Is this reflection a systematic absence?
   *
   * @return a boolean.
   */
  public boolean sysAbs() {
    return (epsilon == 0);
  }

  /**
   * quadForm
   *
   * @param mat an array of double.
   * @return a double.
   */
  public double quadForm(double[][] mat) {
    return h * (h * mat[0][0] + 2 * (k * mat[0][1] + l * mat[0][2]))
        + k * (k * mat[1][1] + 2 * (l * mat[1][2]))
        + l * l * mat[2][2];
  }

  /**
   * getAllowed
   *
   * @return a double.
   */
  public double getAllowed() {
    return ((double) allowed) * (PI / ndiv);
  }

  /**
   * setAllowed
   *
   * @param allowed an int.
   */
  public void setAllowed(int allowed) {
    this.allowed = allowed;
  }

  /**
   * The bin number of this reflection, which is used for resolution dependent R/Rfree.
   *
   * @return the bin number of this reflection.
   */
  public int getBin() {
    return bin;
  }

  /**
   * Set the bin number of this reflection, which is used for resolution dependent R/Rfree.
   *
   * @param bin the bin number of this reflection.
   */
  public void setBin(int bin) {
    this.bin = bin;
  }

  /**
   * Is this reflection centric?
   *
   * @return a boolean.
   */
  public boolean centric() {
    return (allowed != 255);
  }

  /**
   * getEpsilon
   *
   * @return an int.
   */
  public int getEpsilon() {
    return epsilon;
  }

  /**
   * setEpsilon
   *
   * @param eps an int.
   */
  public void setEpsilon(int eps) {
    this.epsilon = eps;
  }

  /**
   * epsilonc
   *
   * @return an int.
   */
  public int epsilonc() {
    if (centric()) {
      return 2 * epsilon;
    } else {
      return epsilon;
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    HKL hkl = (HKL) o;
    return (h == hkl.getH() && k == hkl.getK() && l == hkl.getL());
  }

  /**
   * The h-index of the reflection.
   *
   * @return The h-index of the reflection.
   */
  public int getH() {
    return h;
  }

  /**
   * Set the h-index of the reflection.
   *
   * @param h The h-index of the reflection.
   */
  public void setH(int h) {
    this.h = h;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int hashCode() {
    return Objects.hash(h, k, l);
  }

  /**
   * Get the index of this reflection.
   *
   * @return an int.
   */
  public int getIndex() {
    return index;
  }

  /**
   * Set the index of this reflection.
   *
   * @param index The unique index of this reflection.
   */
  public void setIndex(int index) {
    this.index = index;
  }

  /**
   * Get the k-index of the reflection.
   *
   * @return an int.
   */
  public int getK() {
    return k;
  }

  /**
   * Set the k-index of the reflection.
   *
   * @param k an int.
   */
  public void setK(int k) {
    this.k = k;
  }

  /**
   * Get the l-index of the reflection.
   *
   * @return an int.
   */
  public int getL() {
    return l;
  }

  /**
   * Set the l-index of the reflection.
   *
   * @param l an int.
   */
  public void setL(int l) {
    this.l = l;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String toString() {
    return h + " " + k + " " + l
        + "(allowed: " + allowed + " eps: " + epsilon + ") ";
  }
}
