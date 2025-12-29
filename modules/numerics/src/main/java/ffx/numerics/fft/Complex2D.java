// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import java.util.Random;

import static java.lang.Math.max;
import static java.lang.Math.min;

/**
 * Compute the 2D FFT of complex, double precision input of arbitrary dimensions
 * via 1D Mixed Radix FFTs.
 * <p>
 * For interleaved data, the location of the input point [x, y] within the input array must be: <br>
 * double real = input[x*nextX + y*nextY]<br>
 * double imag = input[x*nextX + y*nextY + im]<br>
 * where <br>
 * nextX = 2
 * nextY = 2*nX
 * im = 1
 * <p>
 * For blocked data along x, the location of the input point [x, y] within the input array must be: <br>
 * double real = input[x*nextX + y*nextY]<br>
 * double imag = input[x*nextX + y*nextY + im]<br>
 * where for BLOCKED_X <br>
 * nextX = 1<br>
 * nextY = 2*nX<br>
 * im = nX<br>
 * and for BLOCKED_XY <br>
 * nextX = 1<br>
 * nextY = nX<br>
 * im = nX*nY<br>
 *
 * @author Michal J. Schnieders
 * @see Complex
 * @since 1.0
 */
public class Complex2D {

  /**
   * The 2D data layout.
   */
  private final DataLayout2D layout;
  /**
   * The external imaginary offset.
   */
  private final int externalIm;
  /**
   * The X-dimension.
   */
  private final int nX;
  /**
   * The Y-dimension.
   */
  private final int nY;
  /**
   * The next real value along the X dimension.
   */
  private final int nextX;
  /**
   * The next real value along the Y dimension.
   */
  private final int nextY;
  /**
   * Compute FFTs along X one at a time.
   */
  private final Complex fftX;
  /**
   * Compute FFTs along Y one at a time.
   */
  private final Complex fftY;
  /**
   * If true, pack FFTs along X (or Y) into a contiguous array to compute all FFTs along X (or Y) at once.
   */
  private boolean packFFTs;
  /**
   * If true, use SIMD instructions.
   */
  private boolean useSIMD;
  /**
   * Number of FFTs to compute along X in one batch. The optimal size depends on cache layout
   * and SIMD register width. The default is to compute 8 at a time (tileSizeX = 8).
   */
  private int nFFTX;
  /**
   * Number of FFTs to compute along Y in one batch. The optimal size depends on cache layout
   * and SIMD register width. The default is to compute 8 at a time (tileSizeY = 8).
   */
  private int nFFTY;
  /**
   * Compute tileSizeX FFTs along the X dimension all at once.
   */
  private final Complex fftXTile;
  /**
   * Compute tileSizeY FFTs along the Y dimension all at once.
   */
  private final Complex fftYTile;
  /**
   * Working array for packed FFTs.
   */
  private final double[] tile;
  /**
   * The offset between real values in the packed data.
   */
  private final int ii;
  /**
   * The offset between any real value and its corresponding imaginary value for the packed data.
   */
  private final int im;
  /**
   * For FFTs along X, the offset between any real value and its corresponding
   * imaginary value for tiled data.
   * This will be 1 for interleaved or nX for blocked.
   */
  private final int imX;
  /**
   * For FFTs along Y, the offset between any real value and its corresponding
   * imaginary value for tiled data.
   * This will be 1 for interleaved or nY for blocked.
   */
  private final int imY;
  /**
   * The offset between real values along the X-dimension in the transposed packed data.
   */
  private final int trNextX;
  /**
   * The offset between real values along the Y-dimension in the transposed packed data.
   */
  private final int trNextY;

  /**
   * Create a new 2D Complex FFT for interleaved data.
   *
   * @param nX The number of points in the X dimension.
   * @param nY The number of points in the Y dimension.
   */
  public Complex2D(int nX, int nY) {
    this(nX, nY, DataLayout2D.INTERLEAVED, 1);
  }

  /**
   * Create a new 2D Complex FFT.
   *
   * @param nX       The number of points in the X dimension.
   * @param nY       The number of points in the Y dimension.
   * @param layout   The data layout.
   * @param imOffset The offset between real and imaginary values.
   */
  public Complex2D(int nX, int nY, DataLayout2D layout, int imOffset) {
    this.nX = nX;
    this.nY = nY;
    this.externalIm = imOffset;
    this.layout = layout;

    int dX = nY;
    String pX = System.getProperty("fft.tileX", Integer.toString(dX));
    try {
      nFFTX = Integer.parseInt(pX);
      if (nFFTX < 1 || nFFTX > nY) {
        nFFTX = dX;
      }
    } catch (Exception e) {
      nFFTX = dX;
    }

    int dY = nX;
    String pY = System.getProperty("fft.tileY", Integer.toString(dY));
    try {
      nFFTY = Integer.parseInt(pY);
      if (nFFTY < 1 || nFFTY > nX) {
        nFFTY = dY;
      }
    } catch (Exception e) {
      nFFTY = dY;
    }

    if (layout == DataLayout2D.INTERLEAVED) {
      im = 1;
      imX = 1;
      imY = 1;
      ii = 2;
      nextX = 2;
      nextY = 2 * nX;
      trNextY = 2;
      trNextX = 2 * nY;
    } else if (layout == DataLayout2D.BLOCKED_X) {
      im = nX;
      if (nFFTX == nX) {
        // All at once.
        imX = nX;
      } else {
        // Internal tiles using BLOCKED_XY
        imX = nX * nFFTX;
      }
      if (nFFTY == nY) {
        // All at once.
        imY = nX;
      } else {
        // Internal Tiles using BLOCKED_XY
        imY = nY * nFFTY;
      }
      ii = 1;
      nextX = 1;
      nextY = 2 * nX;
      trNextY = 1;
      trNextX = 2 * nY;
    } else if (layout == DataLayout2D.BLOCKED_XY) {
      im = nX * nY;
      imX = nX * nFFTX;   // Packed, FFT along X
      imY = nY * nFFTY;   // Packed, FFT along Y
      ii = 1;
      nextX = 1;
      nextY = nX;
      trNextY = 1;
      trNextX = nY;
    } else {
      throw new IllegalArgumentException(" Unsupported data layout: " + layout);
    }

    if (this.externalIm != im) {
      throw new IllegalArgumentException(" Unsupported im offset: " + imOffset);
    }

    // Use SIMD by default.
    useSIMD = true;
    String simd = System.getProperty("fft.simd", Boolean.toString(useSIMD));
    try {
      useSIMD = Boolean.parseBoolean(simd);
    } catch (Exception e) {
      useSIMD = false;
    }

    packFFTs = true;
    String pack = System.getProperty("fft.pack", Boolean.toString(packFFTs));
    try {
      packFFTs = Boolean.parseBoolean(pack);
    } catch (Exception e) {
      packFFTs = false;
    }

    DataLayout1D layout1D;
    if (layout == DataLayout2D.INTERLEAVED) {
      layout1D = DataLayout1D.INTERLEAVED;
    } else {
      layout1D = DataLayout1D.BLOCKED;
    }

    // Create 1D FFTs that will be used for computing 1 FFT at a time.
    fftX = new Complex(nX, layout1D, im);
    fftX.setUseSIMD(useSIMD);
    fftY = new Complex(nY, layout1D, im);
    fftY.setUseSIMD(useSIMD);

    // Create 1D FFTs that will be used for computing more than 1 FFT at a time.
    fftXTile = new Complex(nX, layout1D, imX, nFFTX);
    fftXTile.setUseSIMD(useSIMD);
    fftYTile = new Complex(nY, layout1D, imY, nFFTY);
    fftYTile.setUseSIMD(useSIMD);

    int arraySize = max(nX * nFFTX, nFFTY * nY);
    tile = new double[2 * arraySize];
  }

  /**
   * Get the Data Layout.
   *
   * @return The layout.
   */
  public DataLayout2D getLayout() {
    return layout;
  }

  /**
   * Set the 2D transform to use SIMD instructions.
   *
   * @param useSIMD True to use SIMD instructions.
   */
  public void setUseSIMD(boolean useSIMD) {
    this.useSIMD = useSIMD;
    fftX.setUseSIMD(useSIMD);
    fftY.setUseSIMD(useSIMD);
    fftXTile.setUseSIMD(useSIMD);
    fftYTile.setUseSIMD(useSIMD);
  }

  /**
   * Set the 2D transform to pack FFTs into a contiguous array to compute all FFTs at once.
   *
   * @param packFFTs True to pack FFTs.
   */
  public void setPackFFTs(boolean packFFTs) {
    this.packFFTs = packFFTs;
  }

  /**
   * Compute the 2D FFT.
   *
   * @param input The input array must be of size 2 * nX * nY.
   * @param index The offset into the input array of the first element.
   */
  public void fft(final double[] input, int index) {
    if (!packFFTs) {
      for (int offset = index, y = 0; y < nY; y++, offset += nextY) {
        fftX.fft(input, offset, nextX);
      }
      for (int offset = index, x = 0; x < nX; x++, offset += nextX) {
        fftY.fft(input, offset, nextY);
      }
    } else {
      if (nFFTY == nX) {
        fftYTile.fft(input, index, ii);
      } else {
        for (int x = 0; x < nX; x += nFFTY) {
          getFFTYTile(input, index, tile, x, nFFTY);
          fftYTile.fft(tile, 0, ii);
          setFFTYTile(input, index, tile, x, nFFTY);
        }
      }
      if (nFFTX == nY) {
        transpose(input, index);
        fftXTile.fft(tile, 0, ii);
        unTranspose(input, index);
      } else {
        for (int y = 0; y < nY; y += nFFTX) {
          getFFTXTile(input, index, tile, y, nFFTX);
          fftXTile.fft(tile, 0, ii);
          setFFTXTile(input, index, tile, y, nFFTX);
        }
      }
    }
  }

  /**
   * Compute the 2D IFFT.
   *
   * @param input The input array must be of size 2 * nX * nY.
   * @param index The offset into the input array of the first element.
   */
  public void ifft(final double[] input, int index) {
    if (!packFFTs) {
      for (int offset = index, y = 0; y < nY; y++, offset += nextY) {
        fftX.ifft(input, offset, nextX);
      }
      for (int offset = index, x = 0; x < nX; x++, offset += nextX) {
        fftY.ifft(input, offset, nextY);
      }
    } else {
      if (nFFTY == nX) {
        fftYTile.ifft(input, index, ii);
      } else {
        for (int x = 0; x < nX; x += nFFTY) {
          getFFTYTile(input, index, tile, x, nFFTY);
          fftYTile.ifft(tile, 0, ii);
          setFFTYTile(input, index, tile, x, nFFTY);
        }
      }
      if (nFFTX == nY) {
        transpose(input, index);
        fftXTile.ifft(tile, 0, ii);
        unTranspose(input, index);
      } else {
        for (int y = 0; y < nY; y += nFFTX) {
          getFFTXTile(input, index, tile, y, nFFTX);
          fftXTile.ifft(tile, 0, ii);
          setFFTXTile(input, index, tile, y, nFFTX);
        }
      }
    }
  }

  /**
   * The input is of size nX x nY.
   * Extract a tile of size tileSize x nY.
   * <p>
   * Input order:
   * real(x,y) = input[offset + x*nextX + y*nextY]
   * imag(x,y) = input[offset + x*nextX + y*nextY + externalIm]
   * Output::
   * For each Y, all X-values for the tile are contiguous in memory.
   *
   * @param input    The input data.
   * @param offset   The offset into the input data.
   * @param tile     The output tile.
   * @param firstX   The first x-index of the tile.
   * @param tileSize The size of the tile along X.
   */
  private void getFFTYTile(final double[] input, int offset, final double[] tile, int firstX, int tileSize) {
    int index = 0;
    // Handle padding
    int padX = firstX + tileSize;
    int lastX = min(nX, padX);
    // Outer loop over the Y dimension.
    for (int y = 0; y < nY; y++) {
      int dy = offset + y * nextY;
      // Inner loop over the X dimension.
      for (int x = firstX; x < lastX; x++) {
        int dx = x * nextX;
        double real = input[dx + dy];
        double imag = input[dx + dy + externalIm];
        // Contiguous storage into a packed array.
        tile[index] = real;
        tile[index + imY] = imag;
        index += ii;
      }
      for (int x = lastX; x < padX; x++) {
        index += ii;
      }
    }
  }

  /**
   * Set a tile of size tileSize x nY.
   * Put it into the input array of size nX x nY.
   * <p>
   * Input array order:
   * real(x,y) = input[offset + x*nextX + y*nextY]
   * imag(x,y) = input[offset + x*nextX + y*nextY + externalIm]
   *
   * @param input    The input data.
   * @param offset   The offset into the input data.
   * @param tile     The output tile.
   * @param firstX   The first x-index of the tile.
   * @param tileSize The size of the tile along X.
   */
  private void setFFTYTile(final double[] input, int offset, final double[] tile, int firstX, int tileSize) {
    int index = 0;
    // Handle padding
    int padX = firstX + tileSize;
    int lastX = min(nX, padX);
    // Outer loop over the Y dimension.
    for (int y = 0; y < nY; y++) {
      int dy = offset + y * nextY;
      // Inner loop over the X dimension.
      for (int x = firstX; x < lastX; x++) {
        double real = tile[index];
        double imag = tile[index + imY];
        index += ii;
        int dx = x * nextX;
        input[dx + dy] = real;
        input[dx + dy + externalIm] = imag;
      }
      for (int x = lastX; x < padX; x++) {
        index += ii;
      }
    }
  }

  /**
   * The input is of size nX x nY.
   * Extract a tile of size nX * tileSize.
   * <p>
   * Input order:
   * real(x,y) = input[offset + x*nextX + y*nextY]
   * imag(x,y) = input[offset + x*nextX + y*nextY + externalIm]
   * Tile order:
   * For each X, all Y-values for the tile a contiguous in memory.
   *
   * @param input    The input data.
   * @param offset   The offset into the input data.
   * @param tile     The output tile.
   * @param firstY   The first y-index of the tile.
   * @param tileSize The size of the tile along Y.
   */
  private void getFFTXTile(final double[] input, int offset, final double[] tile, int firstY, int tileSize) {
    int index = 0;
    // Handle padding
    int padY = firstY + tileSize;
    int lastY = min(nY, padY);
    // Outer loop over the X dimension.
    for (int x = 0; x < nX; x++) {
      int dx = offset + x * nextX;
      // Inner loop over the Y dimension.
      for (int y = firstY; y < lastY; y++) {
        int dy = y * nextY;
        double real = input[dx + dy];
        double imag = input[dx + dy + externalIm];
        // Contiguous storage into a packed array.
        tile[index] = real;
        tile[index + imX] = imag;
        index += ii;
      }
      for (int y = lastY; y < padY; y++) {
        index += ii;
      }
    }
  }

  /**
   * Set a tile of size nX * tileSize.
   * Results are written back into the "input" array of size nX x nY.
   * <p>
   * Tile order:
   * For each X, all Y-values for the tile a contiguous in memory.
   * <p>
   * Input order:
   * real(x,y) = input[offset + x*nextX + y*nextY]
   * imag(x,y) = input[offset + x*nextX + y*nextY + externalIm]
   *
   * @param input    The input data.
   * @param offset   The offset into the input data.
   * @param tile     The output tile.
   * @param firstY   The first y-index of the tile.
   * @param tileSize The size of the tile along Y.
   */
  private void setFFTXTile(final double[] input, int offset, final double[] tile, int firstY, int tileSize) {
    int index = 0;
    // Handle padding
    int padY = firstY + tileSize;
    int lastY = min(nY, padY);
    // Outer loop over the X dimension.
    for (int x = 0; x < nX; x++) {
      int dx = offset + x * nextX;
      // Inner loop over the Y dimension.
      for (int y = firstY; y < lastY; y++) {
        double real = tile[index];
        double imag = tile[index + imX];
        index += ii;
        int dy = y * nextY;
        input[dx + dy] = real;
        input[dx + dy + externalIm] = imag;
      }
      for (int y = lastY; y < padY; y++) {
        index += ii;
      }
    }
  }

  /**
   * Pack the input array for Fourier transforms along the X dimension.
   * <p>
   * Input order:
   * real(x,y) = input[offset + x*nextX + y*nextY]
   * imag(x,y) = input[offset + x*nextX + y*nextY + im]
   * Output order:
   * real(x,y) = packedData[y*trNextY + x*trNextX]
   * imag(x,y) = packedData[y*trNextY + x*trXextX + im]
   *
   * @param input  The input data.
   * @param offset The offset into the input data.
   */
  private void transpose(final double[] input, int offset) {
    int index = 0;
    // Outer loop over the X dimension.
    for (int x = 0; x < nX; x++) {
      int dx = offset + x * nextX;
      // Inner loop over the Y dimension (the number of FFTs).
      for (int y = 0; y < nY; y++) {
        double real = input[dx + y * nextY];
        double imag = input[dx + y * nextY + externalIm];
        // Contiguous storage into the packed array.
        tile[index] = real;
        tile[index + im] = imag;
        index += ii;
      }
    }
  }

  /**
   * Unpack the output array after Fourier transforms.
   * <p>
   * Input order:
   * real_xy = packedData[y*trNextY + x*trNextX]
   * imag_xy = packedData[y*trNextY + x*trXextX + im]
   * Output order:
   * real_xy = output[offset + x*nextX + y*nextY]
   * imag_xy = output[offset + x*nextX + y*nextY + im]
   *
   * @param output The output data.
   * @param offset The offset into the output data.
   */
  private void unTranspose(final double[] output, int offset) {
    int index = offset;
    // Outer loop over the Y dimension.
    for (int y = 0; y < nY; y++) {
      int dy = y * trNextY;
      // Inner loop over the X dimension.
      for (int x = 0; x < nX; x++) {
        double real = tile[dy + x * trNextX];
        double imag = tile[dy + x * trNextX + im];
        // Contiguous storage into the output array.
        output[index] = real;
        output[index + externalIm] = imag;
        index += ii;
      }
    }
  }

  /**
   * Test the Complex2D FFT.
   *
   * @param args an array of {@link java.lang.String} objects.
   * @throws java.lang.Exception if any.
   * @since 1.0
   */
  public static void main(String[] args) throws Exception {
    int dimNotFinal = 128;
    int reps = 5;
    boolean blocked = false;
    try {
      dimNotFinal = Integer.parseInt(args[0]);
      if (dimNotFinal < 1) {
        dimNotFinal = 128;
      }
      reps = Integer.parseInt(args[1]);
      if (reps < 1) {
        reps = 5;
      }
      blocked = Boolean.parseBoolean(args[2]);
    } catch (Exception e) {
      //
    }
    final int dim = dimNotFinal;
    System.out.printf("Initializing a %d cubed grid.\n"
            + "The best timing out of %d repetitions will be used.%n",
        dim, reps);
    // One dimension of the serial array divided by the number of threads.
    Complex2D complex2D;
    if (blocked) {
      complex2D = new Complex2D(dim, dim, DataLayout2D.BLOCKED_X, dim);
    } else {
      complex2D = new Complex2D(dim, dim);
    }
    ParallelTeam parallelTeam = new ParallelTeam();
    final double[] data = initRandomData(dim, parallelTeam);

    double toSeconds = 0.000000001;
    long forwardTime = Long.MAX_VALUE;
    long inverseTime = Long.MAX_VALUE;

    // Warm-up
    System.out.println(" Warm Up FFT");
    complex2D.fft(data, 0);
    System.out.println(" Warm Up IFFT");
    complex2D.ifft(data, 0);

    for (int i = 0; i < reps; i++) {
      System.out.printf(" Iteration %d%n", i + 1);
      long time = System.nanoTime();
      complex2D.fft(data, 0);
      time = (System.nanoTime() - time);
      System.out.printf("  FFT:   %9.6f (sec)%n", toSeconds * time);
      if (time < forwardTime) {
        forwardTime = time;
      }
      time = System.nanoTime();
      complex2D.ifft(data, 0);
      time = (System.nanoTime() - time);
      System.out.printf("  IFFT:  %9.6f (sec)%n", toSeconds * time);
      if (time < inverseTime) {
        inverseTime = time;
      }
    }

    System.out.printf(" Best FFT Time:    %9.6f (sec)%n", toSeconds * forwardTime);
    System.out.printf(" Best IFFT Time:   %9.6f (sec)%n", toSeconds * inverseTime);
    parallelTeam.shutdown();
  }

  /**
   * Initialize a 2D data for testing purposes.
   *
   * @param dim The dimension of the 2D grid.
   * @return The 3D data.
   * @since 1.0
   */
  public static double[] initRandomData(int dim, ParallelTeam parallelTeam) {
    int n = dim * dim;
    double[] data = new double[2 * n];
    try {
      parallelTeam.execute(
          new ParallelRegion() {
            @Override
            public void run() {
              try {
                execute(
                    0,
                    dim - 1,
                    new IntegerForLoop() {
                      @Override
                      public void run(final int lb, final int ub) {
                        Random randomNumberGenerator = new Random(1);
                        int index = dim * lb * 2;
                        for (int i = lb; i <= ub; i++) {
                          for (int j = 0; j < dim; j++) {
                            double randomNumber = randomNumberGenerator.nextDouble();
                            data[index] = randomNumber;
                            index += 2;
                          }
                        }
                      }
                    });
              } catch (Exception e) {
                System.out.println(e.getMessage());
                System.exit(-1);
              }
            }
          });
    } catch (Exception e) {
      System.out.println(e.getMessage());
      System.exit(-1);
    }
    return data;
  }
}
