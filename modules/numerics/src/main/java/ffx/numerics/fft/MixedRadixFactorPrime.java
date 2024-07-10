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

import static java.lang.Math.fma;

public class MixedRadixFactorPrime extends MixedRadixFactor {

  public MixedRadixFactorPrime(PassConstants passConstants) {
    super(passConstants);
  }

  /**
   * Handle a general prime number.
   */
  protected void passScalar() {
    if (im != 1) {
      throw new IllegalArgumentException(" Support for large prime factors requires interleaved data.");
    }
    final int dataOffset = i;
    final int retOffset = j;
    final int jump = (factor - 1) * innerLoopLimit;
    for (int i = 0; i < nextInput; i++) {
      ret[retOffset + 2 * i] = data[dataOffset + 2 * i];
      ret[retOffset + 2 * i + 1] = data[dataOffset + 2 * i + 1];
    }
    for (int e = 1; e < (factor - 1) / 2 + 1; e++) {
      for (int i = 0; i < nextInput; i++) {
        int idx = i + e * nextInput;
        int idxc = i + (factor - e) * nextInput;
        ret[retOffset + 2 * idx] = data[dataOffset + 2 * idx] + data[dataOffset + 2 * idxc];
        ret[retOffset + 2 * idx + 1] = data[dataOffset + 2 * idx + 1] + data[dataOffset + 2 * idxc + 1];
        ret[retOffset + 2 * idxc] = data[dataOffset + 2 * idx] - data[dataOffset + 2 * idxc];
        ret[retOffset + 2 * idxc + 1] = data[dataOffset + 2 * idx + 1] - data[dataOffset + 2 * idxc + 1];
      }
    }
    for (int i = 0; i < nextInput; i++) {
      data[dataOffset + 2 * i] = ret[retOffset + 2 * i];
      data[dataOffset + 2 * i + 1] = ret[retOffset + 2 * i + 1];
    }
    for (int e1 = 1; e1 < (factor - 1) / 2 + 1; e1++) {
      for (int i = 0; i < nextInput; i++) {
        int i1 = retOffset + 2 * (i + e1 * nextInput);
        data[dataOffset + 2 * i] += ret[i1];
        data[dataOffset + 2 * i + 1] += ret[i1 + 1];
      }
    }
    double[] twiddl = twiddles[outerLoopLimit];
    for (int e = 1; e < (factor - 1) / 2 + 1; e++) {
      int idx = e;
      double wr, wi;
      int em = e * nextInput;
      int ecm = (factor - e) * nextInput;
      for (int i = 0; i < nextInput; i++) {
        data[dataOffset + 2 * (i + em)] = ret[retOffset + 2 * i];
        data[dataOffset + 2 * (i + em) + 1] = ret[retOffset + 2 * i + 1];
        data[dataOffset + 2 * (i + ecm)] = ret[retOffset + 2 * i];
        data[dataOffset + 2 * (i + ecm) + 1] = ret[retOffset + 2 * i + 1];
      }
      for (int e1 = 1; e1 < (factor - 1) / 2 + 1; e1++) {
        if (idx == 0) {
          wr = 1;
          wi = 0;
        } else {
          wr = twiddl[2 * (idx - 1)];
          wi = -sign * twiddl[2 * (idx - 1) + 1];
        }
        for (int i = 0; i < nextInput; i++) {
          int i1 = retOffset + 2 * (i + e1 * nextInput);
          int i2 = retOffset + 2 * (i + (factor - e1) * nextInput);
          double ap = wr * ret[i1];
          double am = wi * ret[i2 + 1];
          double bp = wr * ret[i1 + 1];
          double bm = wi * ret[i2];
          data[dataOffset + 2 * (i + em)] += (ap - am);
          data[dataOffset + 2 * (i + em) + 1] += (bp + bm);
          data[dataOffset + 2 * (i + ecm)] += (ap + am);
          data[dataOffset + 2 * (i + ecm) + 1] += (bp - bm);
        }
        idx += e;
        idx %= factor;
      }
    }
    /* k = 0 */
    for (int k1 = 0; k1 < innerLoopLimit; k1++) {
      ret[retOffset + 2 * k1] = data[dataOffset + 2 * k1];
      ret[retOffset + 2 * k1 + 1] = data[dataOffset + 2 * k1 + 1];
    }
    for (int e1 = 1; e1 < factor; e1++) {
      for (int k1 = 0; k1 < innerLoopLimit; k1++) {
        int i = retOffset + 2 * (k1 + e1 * innerLoopLimit);
        int i1 = dataOffset + 2 * (k1 + e1 * nextInput);
        ret[i] = data[i1];
        ret[i + 1] = data[i1 + 1];
      }
    }
    int i = innerLoopLimit;
    int j = product;
    for (int k = 1; k < outerLoopLimit; k++) {
      for (int k1 = 0; k1 < innerLoopLimit; k1++) {
        ret[retOffset + 2 * j] = data[dataOffset + 2 * i];
        ret[retOffset + 2 * j + 1] = data[dataOffset + 2 * i + 1];
        i++;
        j++;
      }
      j += jump;
    }
    i = innerLoopLimit;
    j = product;
    for (int k = 1; k < outerLoopLimit; k++) {
      twiddl = twiddles[k];
      for (int k1 = 0; k1 < innerLoopLimit; k1++) {
        for (int e1 = 1; e1 < factor; e1++) {
          int i1 = dataOffset + 2 * (i + e1 * nextInput);
          double xr = data[i1];
          double xi = data[i1 + 1];
          double wr = twiddl[2 * (e1 - 1)];
          double wi = -sign * twiddl[2 * (e1 - 1) + 1];
          int i2 = retOffset + 2 * (j + e1 * innerLoopLimit);
          ret[i2] = fma(wr, xr, -wi * xi);
          ret[i2 + 1] = fma(wr, xi, wi * xr);
        }
        i++;
        j++;
      }
      j += jump;
    }
  }

  /**
   * There is no SIMD support for general prime numbers, so the scalar method is used.
   */
  protected void passSIMD() {
    passScalar();
  }

}
