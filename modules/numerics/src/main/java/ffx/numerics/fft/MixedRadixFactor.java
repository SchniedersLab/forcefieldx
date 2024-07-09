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
package ffx.numerics.fft;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorShuffle;
import jdk.incubator.vector.VectorSpecies;

public abstract class MixedRadixFactor {

  /**
   * The preferred SIMD vector size.
   */
  protected static final VectorSpecies<Double> DOUBLE_SPECIES = DoubleVector.SPECIES_PREFERRED;
  /**
   * Vector used to change the sign of the imaginary members of the vector via multiplication.
   */
  protected static final DoubleVector negateIm;
  /**
   * Vector used to change the sign of the real members of the vector via multiplication.
   */
  protected static final DoubleVector negateRe;
  /**
   * Shuffle used to swap real and imaginary members of the vector.
   */
  protected static final VectorShuffle<Double> shuffleReIm;
  /**
   * The number of contiguous elements that will be read from the input data array.
   */
  protected static final int SPECIES_LENGTH = DOUBLE_SPECIES.length();
  /**
   * The number of complex elements that will be processed in each inner loop iteration.
   * The number of elements to process in the inner loop must be evenly divisible by this loop increment.
   */
  protected static final int LOOP_INCREMENT = SPECIES_LENGTH / 2;

  static {
    // Assume that 512 is the largest vector size.
    if (SPECIES_LENGTH > 8) {
      throw new IllegalStateException("Unsupported SIMD vector size: " + SPECIES_LENGTH);
    }
    double[] negateReal = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    int[] shuffleMask = {1, 0, 3, 2, 5, 4, 7, 6};
    negateRe = DoubleVector.fromArray(DOUBLE_SPECIES, negateReal, 0);
    negateIm = negateRe.mul(-1.0);
    shuffleReIm = VectorShuffle.fromArray(DOUBLE_SPECIES, shuffleMask, 0);
  }

  protected final int factor;
  protected final int product;
  protected final int outerLoopLimit;
  protected final int innerLoopLimit;
  protected final int nextInput;
  protected final int di;
  protected final int dj;
  protected final double[][] twiddles;
  protected final int jstep;

  public MixedRadixFactor(PassConstants passConstants) {
    this.factor = passConstants.factor();
    this.product = passConstants.product();
    this.outerLoopLimit = passConstants.outerLoopLimit();
    this.innerLoopLimit = passConstants.innerLoopLimit();
    this.nextInput = passConstants.nextInput();
    this.di = passConstants.di();
    this.dj = passConstants.dj();
    this.twiddles = passConstants.twiddles();
    this.jstep = (factor - 1) * dj;
  }

  protected abstract void pass(PassData passData);

  protected abstract void passSIMD(PassData passData);
}
