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

import static jdk.incubator.vector.DoubleVector.SPECIES_128;
import static jdk.incubator.vector.DoubleVector.SPECIES_256;
import static jdk.incubator.vector.DoubleVector.SPECIES_512;

public class MixedRadixFactor2 extends MixedRadixFactor {

  public MixedRadixFactor2(PassConstants passConstants) {
    super(passConstants);
  }

  /**
   * Handle factors of 2.
   */
  @Override
  protected void passScalar() {
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
      final double z0_r = data[i];
      final double z0_i = data[i + im];
      final int idi = i + di;
      final double z1_r = data[idi];
      final double z1_i = data[idi + im];
      ret[j] = z0_r + z1_r;
      ret[j + im] = z0_i + z1_i;
      final double x_r = z0_r - z1_r;
      final double x_i = z0_i - z1_i;
      final int jdj = j + dj;
      ret[jdj] = x_r;
      ret[jdj + im] = x_i;
    }
    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      final double w_r = twids[0];
      final double w_i = -sign * twids[1];
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
        final double z0_r = data[i];
        final double z0_i = data[i + im];
        final int idi = i + di;
        final double z1_r = data[idi];
        final double z1_i = data[idi + im];
        ret[j] = z0_r + z1_r;
        ret[j + im] = z0_i + z1_i;
        final int jdj = j + dj;
        multiplyAndStore(z0_r - z1_r, z0_i - z1_i, w_r, w_i, ret, jdj, jdj + im);
      }
    }
  }

  /**
   * Handle factors of 2 using SIMD vectors.
   */
  @Override
  protected void passSIMD() {
    if (im == 1) {
      interleaved();
    } else {
      blocked();
    }
  }

  private void interleaved() {
    if (innerLoopLimit % LOOP == 0) {
      // Use the preferred SIMD vector.
      switch (LENGTH) {
        case 2:
          interleaved128();
          break;
        case 4:
          interleaved256();
          break;
        case 8:
          interleaved512();
          break;
      }
    } else {
      // If the inner loop limit is not divisible by the loop increment, use largest SIMD vector that fits.
      switch (innerLoopLimit) {
        case 1:
          // 1 Complex
          interleaved128();
          break;
        case 2:
          // 2 Complex
          interleaved256();
          break;
        case 4:
          // 4 Complex
          interleaved512();
          break;
        default:
          throw new IllegalStateException(" Unsupported inner loop limit: " + innerLoopLimit);
      }
    }
  }

  private void blocked() {
    if (innerLoopLimit % BLOCK_LOOP == 0) {
      // The preferred SIMD vector length is a multiple of the inner loop limit and can be used.
      switch (LENGTH) {
        case 2:
          blocked128();
          break;
        case 4:
          blocked256();
          break;
        case 8:
          blocked512();
          break;
      }
    } else {
      // If the inner loop limit is not divisible by the preferred loop increment, use largest SIMD vector that fits.
      switch (innerLoopLimit) {
        case 1:
          // Use the scalar method.
          passScalar();
          break;
        case 2:
          // 2 Real and 2 Imaginary per loop iteration.
          blocked128();
          break;
        case 4:
          // 4 Real and 4 Imaginary per loop iteration.
          blocked256();
          break;
        case 8:
          // 8 Real and 8 Imaginary per loop iteration.
          blocked512();
          break;
        default:
          throw new IllegalStateException(" Unsupported inner loop limit: " + innerLoopLimit);
      }
    }
  }

  /**
   * Handle factors of 2 using the 128-bit SIMD vectors.
   */
  private void blocked128() {
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_128, i += LENGTH_128, j += LENGTH_128) {
      DoubleVector z0_r = DoubleVector.fromArray(SPECIES_128, data, i);
      DoubleVector z1_r = DoubleVector.fromArray(SPECIES_128, data, i + di);
      DoubleVector z0_i = DoubleVector.fromArray(SPECIES_128, data, i + im);
      DoubleVector z1_i = DoubleVector.fromArray(SPECIES_128, data, i + di + im);
      z0_r.add(z1_r).intoArray(ret, j);
      z0_i.add(z1_i).intoArray(ret, j + im);
      z0_r.sub(z1_r).intoArray(ret, j + dj);
      z0_i.sub(z1_i).intoArray(ret, j + dj + im);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector
          w_r = DoubleVector.broadcast(SPECIES_128, twids[0]),
          w_i = DoubleVector.broadcast(SPECIES_128, -sign * twids[1]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_128, i += LENGTH_128, j += LENGTH_128) {
        DoubleVector z0_r = DoubleVector.fromArray(SPECIES_128, data, i);
        DoubleVector z1_r = DoubleVector.fromArray(SPECIES_128, data, i + di);
        DoubleVector z0_i = DoubleVector.fromArray(SPECIES_128, data, i + im);
        DoubleVector z1_i = DoubleVector.fromArray(SPECIES_128, data, i + di + im);
        z0_r.add(z1_r).intoArray(ret, j);
        z0_i.add(z1_i).intoArray(ret, j + im);
        multiplyAndStoreBlocked(z0_r.sub(z1_r), z0_i.sub(z1_i), w_r, w_i, ret, j + dj, j + dj + im);
      }
    }
  }

  /**
   * Handle factors of 2 using the 256-bit SIMD vectors.
   */
  private void blocked256() {
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_256, i += LENGTH_256, j += LENGTH_256) {
      DoubleVector z0_r = DoubleVector.fromArray(SPECIES_256, data, i);
      DoubleVector z1_r = DoubleVector.fromArray(SPECIES_256, data, i + di);
      DoubleVector z0_i = DoubleVector.fromArray(SPECIES_256, data, i + im);
      DoubleVector z1_i = DoubleVector.fromArray(SPECIES_256, data, i + di + im);
      z0_r.add(z1_r).intoArray(ret, j);
      z0_i.add(z1_i).intoArray(ret, j + im);
      z0_r.sub(z1_r).intoArray(ret, j + dj);
      z0_i.sub(z1_i).intoArray(ret, j + dj + im);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector
          w_r = DoubleVector.broadcast(SPECIES_256, twids[0]),
          w_i = DoubleVector.broadcast(SPECIES_256, -sign * twids[1]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_256, i += LENGTH_256, j += LENGTH_256) {
        DoubleVector z0_r = DoubleVector.fromArray(SPECIES_256, data, i);
        DoubleVector z1_r = DoubleVector.fromArray(SPECIES_256, data, i + di);
        DoubleVector z0_i = DoubleVector.fromArray(SPECIES_256, data, i + im);
        DoubleVector z1_i = DoubleVector.fromArray(SPECIES_256, data, i + di + im);
        z0_r.add(z1_r).intoArray(ret, j);
        z0_i.add(z1_i).intoArray(ret, j + im);
        multiplyAndStoreBlocked(z0_r.sub(z1_r), z0_i.sub(z1_i), w_r, w_i, ret, j + dj, j + dj + im);
      }
    }
  }

  /**
   * Handle factors of 2 using the 512-bit SIMD vectors.
   */
  private void blocked512() {
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_512, i += LENGTH_512, j += LENGTH_512) {
      DoubleVector z0_r = DoubleVector.fromArray(SPECIES_512, data, i);
      DoubleVector z1_r = DoubleVector.fromArray(SPECIES_512, data, i + di);
      DoubleVector z0_i = DoubleVector.fromArray(SPECIES_512, data, i + im);
      DoubleVector z1_i = DoubleVector.fromArray(SPECIES_512, data, i + di + im);
      z0_r.add(z1_r).intoArray(ret, j);
      z0_i.add(z1_i).intoArray(ret, j + im);
      z0_r.sub(z1_r).intoArray(ret, j + dj);
      z0_i.sub(z1_i).intoArray(ret, j + dj + im);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector
          w_r = DoubleVector.broadcast(SPECIES_512, twids[0]),
          w_i = DoubleVector.broadcast(SPECIES_512, -sign * twids[1]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_512, i += LENGTH_512, j += LENGTH_512) {
        DoubleVector z0_r = DoubleVector.fromArray(SPECIES_512, data, i);
        DoubleVector z1_r = DoubleVector.fromArray(SPECIES_512, data, i + di);
        DoubleVector z0_i = DoubleVector.fromArray(SPECIES_512, data, i + im);
        DoubleVector z1_i = DoubleVector.fromArray(SPECIES_512, data, i + di + im);
        z0_r.add(z1_r).intoArray(ret, j);
        z0_i.add(z1_i).intoArray(ret, j + im);
        multiplyAndStoreBlocked(z0_r.sub(z1_r), z0_i.sub(z1_i), w_r, w_i, ret, j + dj, j + dj + im);
      }
    }
  }

  /**
   * Handle factors of 2 using the 128-bit SIMD vectors.
   */
  private void interleaved128() {
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_128, i += LENGTH_128, j += LENGTH_128) {
      DoubleVector
          z0 = DoubleVector.fromArray(SPECIES_128, data, i),
          z1 = DoubleVector.fromArray(SPECIES_128, data, i + di);
      z0.add(z1).intoArray(ret, j);
      z0.sub(z1).intoArray(ret, j + dj);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector wr = DoubleVector.broadcast(SPECIES_128, twids[0]);
      DoubleVector wi = DoubleVector.broadcast(SPECIES_128, -sign * twids[1]).mul(NEGATE_IM_128);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_128, i += LENGTH_128, j += LENGTH_128) {
        DoubleVector
            z0 = DoubleVector.fromArray(SPECIES_128, data, i),
            z1 = DoubleVector.fromArray(SPECIES_128, data, i + di);
        z0.add(z1).intoArray(ret, j);
        multiplyAndStore128(z0.sub(z1), wr, wi, ret, j + dj);
      }
    }
  }

  /**
   * Handle factors of 2 using the 256-bit SIMD vectors.
   */
  private void interleaved256() {
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_256, i += LENGTH_256, j += LENGTH_256) {
      DoubleVector
          z0 = DoubleVector.fromArray(SPECIES_256, data, i),
          z1 = DoubleVector.fromArray(SPECIES_256, data, i + di);
      z0.add(z1).intoArray(ret, j);
      z0.sub(z1).intoArray(ret, j + dj);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector wr = DoubleVector.broadcast(SPECIES_256, twids[0]);
      DoubleVector wi = DoubleVector.broadcast(SPECIES_256, -sign * twids[1]).mul(NEGATE_IM_256);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_256, i += LENGTH_256, j += LENGTH_256) {
        DoubleVector
            z0 = DoubleVector.fromArray(SPECIES_256, data, i),
            z1 = DoubleVector.fromArray(SPECIES_256, data, i + di);
        z0.add(z1).intoArray(ret, j);
        multiplyAndStore256(z0.sub(z1), wr, wi, ret, j + dj);
      }
    }
  }

  /**
   * Handle factors of 2 using the 512-bit SIMD vectors.
   */
  private void interleaved512() {
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_512, i += LENGTH_512, j += LENGTH_512) {
      DoubleVector
          z0 = DoubleVector.fromArray(SPECIES_512, data, i),
          z1 = DoubleVector.fromArray(SPECIES_512, data, i + di);
      z0.add(z1).intoArray(ret, j);
      z0.sub(z1).intoArray(ret, j + dj);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector wr = DoubleVector.broadcast(SPECIES_512, twids[0]);
      DoubleVector wi = DoubleVector.broadcast(SPECIES_512, -sign * twids[1]).mul(NEGATE_IM_512);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_512, i += LENGTH_512, j += LENGTH_512) {
        DoubleVector
            z0 = DoubleVector.fromArray(SPECIES_512, data, i),
            z1 = DoubleVector.fromArray(SPECIES_512, data, i + di);
        z0.add(z1).intoArray(ret, j);
        multiplyAndStore512(z0.sub(z1), wr, wi, ret, j + dj);
      }
    }
  }

}