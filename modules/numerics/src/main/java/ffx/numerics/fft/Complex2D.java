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

/**
 * Compute the 2D FFT of complex, double precision input of arbitrary dimensions
 * via 1D Mixed Radix FFTs.
 * <p>
 * For interleaved data, the location of the input point [x, y] within the input array must be: <br>
 * double real = input[x*nextX + y*nextY]<br>
 * double imag = input[x*nextX + y*nextY + 1]<br>
 * where <br>
 * int nextX = 2
 * int nextY = 2*nX
 * <p>
 * For blocked data along x, the location of the input point [x, y] within the input array must be: <br>
 * double real = input[x*nextX + y*nextY]<br>
 * double imag = input[x*nextX + y*nextY + nX]<br>
 * where <br>
 * int nextX = 1
 * int nextY = 2*nX
 *
 * @author Michal J. Schnieders
 * @see Complex
 * @since 1.0
 */
public class Complex2D {

  private final DataLayout1D layout;
  private final int im;
  private final int nX, nY;
  private final int nextX, nextY;
  private final Complex fftX, fftY;

  /**
   * Create a new 2D Complex FFT for interleaved data.
   *
   * @param nX The number of points in the X dimension.
   * @param nY The number of points in the Y dimension.
   */
  public Complex2D(int nX, int nY) {
    this(nX, nY, DataLayout1D.INTERLEAVED, 1);
  }

  /**
   * Create a new 2D Complex FFT.
   *
   * @param nX     The number of points in the X dimension.
   * @param nY     The number of points in the Y dimension.
   * @param layout The data layout.
   */
  public Complex2D(int nX, int nY, DataLayout1D layout, int imOffset) {
    this.nX = nX;
    this.nY = nY;
    this.im = imOffset;
    this.layout = layout;
    if (layout == DataLayout1D.INTERLEAVED) {
      nextX = 2;
      nextY = nextX * nX;
    } else {
      nextX = 1;
      if (im == nX) {
        nextY = 2 * nX;
      } else if (im == nX * nY) {
        nextY = nY * nX;
      } else {
        throw new IllegalArgumentException("Invalid im offset: " + im);
      }
    }
    fftX = new Complex(nX, layout, im);
    if (nX == nY) {
      fftY = fftX;
    } else {
      fftY = new Complex(nY, layout, im);
    }
  }

  /**
   * Get the Data Layout.
   *
   * @return The layout.
   */
  public DataLayout1D getLayout() {
    return layout;
  }

  /**
   * Compute the 2D FFT.
   *
   * @param input The input array must be of size 2 * nX * nY.
   * @param index The offset into the input array of the first element.
   */
  public void fft(final double[] input, int index) {
    for (int offset = index, y = 0; y < nY; y++, offset += nextY) {
      fftX.fft(input, offset, nextX);
    }
    for (int offset = index, x = 0; x < nX; x++, offset += nextX) {
      fftY.fft(input, offset, nextY);
    }
  }

  /**
   * Compute the 2D IFFT.
   *
   * @param input The input array must be of size 2 * nX * nY.
   * @param index The offset into the input array of the first element.
   */
  public void ifft(final double[] input, int index) {
    for (int offset = index, x = 0; x < nX; x++, offset += nextX) {
      fftY.ifft(input, offset, nextY);
    }
    for (int offset = index, y = 0; y < nY; y++, offset += nextY) {
      fftX.ifft(input, offset, nextX);
    }
  }
}
