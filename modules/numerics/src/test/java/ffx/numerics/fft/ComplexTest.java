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

import ffx.utilities.FFXTest;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import java.util.Arrays;
import java.util.Collection;
import java.util.Random;

import static org.junit.Assert.assertEquals;

/**
 * @author Michael J. Schnieders
 */
@RunWith(Parameterized.class)
public class ComplexTest extends FFXTest {

  private final int n;
  private final String info;
  private final boolean preferred;
  private final double[] data;
  private final double[] dataBlocked;
  private final double[] orig;
  private final double[] origBlocked;
  private final double[] dft;

  public ComplexTest(String info, int n, boolean preferred) {
    this.info = info;
    this.n = n;
    this.preferred = preferred;
    data = new double[n * 2];
    dataBlocked = new double[n * 2];
    orig = new double[n * 2];
    origBlocked = new double[n * 2];
    dft = new double[n * 2];
    Random r = new Random();
    for (int i = 0; i < n; i++) {
      double d = r.nextDouble();
      data[i * 2] = d;
      orig[i * 2] = d;
      dataBlocked[i] = d;
      origBlocked[i] = d;
    }
  }

  @Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][]{
            {"Test n = 162 with factors [6, 3, 3, 3]", 162, true},
            {"Test n = 160 with factors [5, 4, 4, 2]", 160, true},
            {"Test n = 120 with factors [6, 5, 4]", 120, true},
            {"Test n = 64 with factors [4, 4, 4]", 64, true},
            {"Test n = 48 with factors [6, 4, 2]", 48, true},
            {"Test n = 21 with factors [7, 3]", 21, true},
            {"Test n = 38 with factors [2, 19]", 38, false},
            {"Test n = 22 with factors [2, 11]", 22, false},
        });
  }

  /**
   * Test of preferredDimension method, of class Complex.
   */
  @Test
  public void testPreferredDimension() {
    boolean result = Complex.preferredDimension(n);
    assertEquals(info, preferred, result);
  }

  /**
   * Test of fft method, of class Complex.
   */
  @Test
  public void testInterleavedScalar() {
    double tolerance = 1.0e-11;

    int offset = 0;
    int stride = 2;
    Complex complex = new Complex(n);
    complex.setUseSIMD(false);

    // System.out.println(info + "\n Factors " + Arrays.toString(complex.getFactors()));

    long dftTime = System.nanoTime();
    Complex.dft(data, dft);
    dftTime = System.nanoTime() - dftTime;
    String dftString = " DFT Time: " + dftTime * 1.0e-9 + " s\n";

    long fftTime = System.nanoTime();
    complex.fft(data, offset, stride);
    fftTime = System.nanoTime() - fftTime;
    String fftString = " FFT Time: " + fftTime * 1.0e-9 + " s";

    // Test the FFT is equals the DFT result.
    for (int i = 0; i < 2 * n; i++) {
      assertEquals(" Forward " + info + " at position: " + i, dft[i], data[i], tolerance);
    }

    // The FFT is faster than the DFT.
    String message = fftString + dftString;
    // assertTrue(message, fftTime < dftTime);

    // Test that X = IFFT(FFT(X)).
    complex.inverse(data, offset, stride);
    for (int i = 0; i < n; i++) {
      double orig = this.orig[i * 2];
      double actual = data[i * 2];
      assertEquals(" IFFT(FFT(X)) " + info + " at position: " + i, orig, actual, tolerance);
    }
  }

  /**
   * Test of fft method, of class Complex.
   */
  @Test
  public void testBlockedScalar() {
    if (!this.preferred) {
      return;
    }
    double tolerance = 1.0e-11;
    int offset = 0;
    int stride = 1;
    Complex complex = new Complex(n, DataLayout1D.BLOCKED, n);
    complex.setUseSIMD(false);

    // System.out.println(info + "\n Factors " + Arrays.toString(complex.getFactors()));

    long dftTime = System.nanoTime();
    Complex.dftBlocked(dataBlocked, dft);
    dftTime = System.nanoTime() - dftTime;
    String dftString = " DFT Time: " + dftTime * 1.0e-9 + " s\n";

    long fftTime = System.nanoTime();
    complex.fft(dataBlocked, offset, stride);
    fftTime = System.nanoTime() - fftTime;
    String fftString = " FFT Time: " + fftTime * 1.0e-9 + " s";

    // Test the FFT is equals the DFT result.
    for (int i = 0; i < 2 * n; i++) {
      assertEquals(" Forward " + info + " at position: " + i, dft[i], dataBlocked[i], tolerance);
    }

    // The FFT is faster than the DFT.
    String message = fftString + dftString;
    // assertTrue(message, fftTime < dftTime);

    // Test that X = IFFT(FFT(X)).
    complex.inverse(dataBlocked, offset, stride);
    for (int i = 0; i < n; i++) {
      double orig = this.origBlocked[i];
      double actual = dataBlocked[i];
      assertEquals(" IFFT(FFT(X)) " + info + " at position: " + i, orig, actual, tolerance);
    }
  }

  /**
   * Test of fft method, of class Complex.
   */
  @Test
  public void testInterleavedSIMD() {
    double tolerance = 1.0e-11;

    int offset = 0;
    int stride = 2;
    Complex complex = new Complex(n);
    complex.setUseSIMD(true);

    // System.out.println(info + "\n Factors " + Arrays.toString(complex.getFactors()));

    long dftTime = System.nanoTime();
    Complex.dft(data, dft);
    dftTime = System.nanoTime() - dftTime;
    String dftString = " DFT Time: " + dftTime * 1.0e-9 + " s\n";

    long fftTime = System.nanoTime();
    complex.fft(data, offset, stride);
    fftTime = System.nanoTime() - fftTime;
    String fftString = " FFT Time: " + fftTime * 1.0e-9 + " s";

    // Test the FFT is equals the DFT result.
    for (int i = 0; i < 2 * n; i++) {
      assertEquals(" Forward " + info + " at position: " + i, dft[i], data[i], tolerance);
    }

    // The FFT is faster than the DFT.
    String message = fftString + dftString;
    // assertTrue(message, fftTime < dftTime);

    // Test that X = IFFT(FFT(X)).
    complex.inverse(data, offset, stride);
    for (int i = 0; i < n; i++) {
      double orig = this.orig[i * 2];
      double actual = data[i * 2];
      assertEquals(" IFFT(FFT(X)) " + info + " at position: " + i, orig, actual, tolerance);
    }
  }

  /**
   * Test of fft method, of class Complex.
   */
  @Test
  public void testBlockedSIMD() {
    if (!this.preferred) {
      return;
    }
    double tolerance = 1.0e-11;
    int offset = 0;
    int stride = 1;
    Complex complex = new Complex(n, DataLayout1D.BLOCKED, n);
    complex.setUseSIMD(true);

    long dftTime = System.nanoTime();
    Complex.dftBlocked(dataBlocked, dft);
    dftTime = System.nanoTime() - dftTime;
    String dftString = " DFT Time: " + dftTime * 1.0e-9 + " s\n";

    long fftTime = System.nanoTime();
    complex.fft(dataBlocked, offset, stride);
    fftTime = System.nanoTime() - fftTime;
    String fftString = " FFT Time: " + fftTime * 1.0e-9 + " s";

    // Test the FFT is equals the DFT result.
    for (int i = 0; i < 2 * n; i++) {
      assertEquals(" Forward " + info + " at position: " + i, dft[i], dataBlocked[i], tolerance);
    }

    // The FFT is faster than the DFT.
    String message = fftString + dftString;
    // assertTrue(message, fftTime < dftTime);

    // Test that X = IFFT(FFT(X)).
    complex.inverse(dataBlocked, offset, stride);
    for (int i = 0; i < n; i++) {
      double orig = this.origBlocked[i];
      double actual = dataBlocked[i];
      assertEquals(" IFFT(FFT(X)) " + info + " at position: " + i, orig, actual, tolerance);
    }
  }

  /**
   * Test of fft method, of class Complex, using N FFTs packed next to each other.
   * <p>
   * This test is designed to mimic a 2D transform, so we'll pack N FFTs into a 1D array of size (2*N)^2
   */
  @Test
  public void testInterleavedScalarNFFT() {
    double tolerance = 1.0e-11;

    int nFFTs = n;
    int nextFFT = 2 * n;
    int offset = 0;
    int stride = 2;
    int im = 1;
    Complex complex = new Complex(n, DataLayout1D.INTERLEAVED, im, nFFTs);
    complex.setUseSIMD(false);

    long dftTime = System.nanoTime();
    Complex.dft(data, dft);
    dftTime = System.nanoTime() - dftTime;
    String dftString = " DFT Time: " + dftTime * 1.0e-9 + " s\n";

    long fftTime = System.nanoTime();
    double[] data2 = new double[2 * n * nFFTs];
    int index = 0;
    for (int ii = 0; ii < 2 * n; ii += stride) {
      for (int i = 0; i < nFFTs; i++) {
        data2[index] = data[ii];
        //data2[index + im] = data[ii + im];
        index += stride;
      }
    }

    complex.fft(data2, offset, stride, nextFFT);
    fftTime = System.nanoTime() - fftTime;
    String fftString = " FFT Time: " + fftTime * 1.0e-9 + " s";

    // Test the FFT is equals the DFT result.
    index = 0;
    for (int ii = 0; ii < n; ii++) {
      for (int i = 0; i < nFFTs; i++) {
        assertEquals(" FFT " + i + " Forward " + info + " at position: " + ii,
            dft[2 * ii], data2[index], tolerance);
        assertEquals(" FFT " + i + " Forward " + info + " at position: " + ii,
            dft[2 * ii + 1], data2[index + 1], tolerance);
        index += 2;
      }
    }

    // The FFT is faster than the DFT.
    String message = fftString + dftString;
    // assertTrue(message, fftTime < dftTime);

    // Test that X = IFFT(FFT(X)).
    complex.inverse(data2, offset, stride, nextFFT);

    // Test the FFT is equals the DFT result.
    index = 0;
    for (int ii = 0; ii < nextFFT; ii += 2) {
      for (int i = 0; i < nFFTs; i++) {
        assertEquals(" IFFT(FFT(X)) " + info + " at position: " + ii + " FFT " + i,
            orig[ii], data2[index], tolerance);
        index += 2;
      }
    }
  }

  /**
   * Test of fft method, of class Complex, using N FFTs packed next to each other.
   * <p>
   * This test is designed to mimic a 2D transform, so we'll pack N FFTs into a 1D array of size (2*N)^2
   */
  @Test
  public void testInterleavedSIMDNFFT() {
    double tolerance = 1.0e-11;

    int nFFTs = n;
    int nextFFT = 2 * n;
    int offset = 0;
    int stride = 2;
    int im = 1;
    Complex complex = new Complex(n, DataLayout1D.INTERLEAVED, im, nFFTs);
    complex.setUseSIMD(true);

    long dftTime = System.nanoTime();
    Complex.dft(data, dft);
    dftTime = System.nanoTime() - dftTime;
    String dftString = " DFT Time: " + dftTime * 1.0e-9 + " s\n";

    long fftTime = System.nanoTime();
    double[] data2 = new double[2 * n * nFFTs];
    int index = 0;
    for (int ii = 0; ii < 2 * n; ii += stride) {
      for (int i = 0; i < nFFTs; i++) {
        data2[index] = data[ii];
        //data2[index + im] = data[ii + im];
        index += stride;
      }
    }

    complex.fft(data2, offset, stride, nextFFT);
    fftTime = System.nanoTime() - fftTime;
    String fftString = " FFT Time: " + fftTime * 1.0e-9 + " s";

    // Test the FFT is equals the DFT result.
    index = 0;
    for (int ii = 0; ii < n; ii++) {
      for (int i = 0; i < nFFTs; i++) {
        assertEquals(" FFT " + i + " Forward " + info + " at position: " + ii,
            dft[2 * ii], data2[index], tolerance);
        assertEquals(" FFT " + i + " Forward " + info + " at position: " + ii,
            dft[2 * ii + 1], data2[index + 1], tolerance);
        index += 2;
      }
    }

    // The FFT is faster than the DFT.
    String message = fftString + dftString;
    // assertTrue(message, fftTime < dftTime);

    // Test that X = IFFT(FFT(X)).
    complex.inverse(data2, offset, stride, nextFFT);

    // Test the FFT is equals the DFT result.
    index = 0;
    for (int ii = 0; ii < nextFFT; ii += 2) {
      for (int i = 0; i < nFFTs; i++) {
        assertEquals(" IFFT(FFT(X)) " + info + " at position: " + ii + " FFT " + i,
            orig[ii], data2[index], tolerance);
        index += 2;
      }
    }
  }

  /**
   * Test of fft method, of class Complex, using N FFTs packed next to each other.
   * <p>
   * This test is designed to mimic a 2D transform, so we'll pack N FFTs into a 1D array of size (2*N)^2
   */
  @Test
  public void testBlockedScalarNFFT() {
    if (!this.preferred) {
      return;
    }
    double tolerance = 1.0e-11;

    int nFFTs = n;
    int nextFFT = 2 * n;
    int offset = 0;
    int stride = 1;
    int im = n * nFFTs;
    Complex complex = new Complex(n, DataLayout1D.BLOCKED, im, nFFTs);
    complex.setUseSIMD(false);

    long dftTime = System.nanoTime();
    Complex.dftBlocked(dataBlocked, dft);
    dftTime = System.nanoTime() - dftTime;
    String dftString = " DFT Time: " + dftTime * 1.0e-9 + " s\n";

    long fftTime = System.nanoTime();
    double[] data2 = new double[2 * n * nFFTs];
    int index = 0;
    for (int ii = 0; ii < n; ii += stride) {
      for (int i = 0; i < nFFTs; i++) {
        data2[index++] = dataBlocked[ii];
      }
    }

    complex.fft(data2, offset, stride, nextFFT);
    fftTime = System.nanoTime() - fftTime;
    String fftString = " FFT Time: " + fftTime * 1.0e-9 + " s";

    // Test the FFT is equals the DFT result.
    index = 0;
    for (int ii = 0; ii < n; ii++) {
      for (int i = 0; i < nFFTs; i++) {
        assertEquals(" FFT " + i + " Forward " + info + " at position: " + ii,
            dft[ii], data2[index], tolerance);
        assertEquals(" FFT " + i + " Forward " + info + " at position: " + ii,
            dft[ii + n], data2[index + im], tolerance);
        index++;
      }
    }

    // The FFT is faster than the DFT.
    String message = fftString + dftString;
    // assertTrue(message, fftTime < dftTime);

    // Test that X = IFFT(FFT(X)).
    complex.inverse(data2, offset, stride, nextFFT);

    // Test the FFT is equals the DFT result.
    index = 0;
    for (int ii = 0; ii < nextFFT; ii += 2) {
      for (int i = 0; i < nFFTs; i++) {
        assertEquals(" IFFT(FFT(X)) " + info + " at position: " + ii + " FFT " + i,
            orig[ii], data2[index], tolerance);
        index++;
      }
    }
  }

  /**
   * Test of fft method, of class Complex, using N FFTs packed next to each other.
   * <p>
   * This test is designed to mimic a 2D transform, so we'll pack N FFTs into a 1D array of size (2*N)^2
   */
  @Test
  public void testBlockedSIMDNFFT() {
    if (!this.preferred) {
      return;
    }
    double tolerance = 1.0e-11;

    int nFFTs = n;
    int nextFFT = 2 * n;
    int offset = 0;
    int stride = 1;
    int im = n * nFFTs;
    Complex complex = new Complex(n, DataLayout1D.BLOCKED, im, nFFTs);
    complex.setUseSIMD(true);

    long dftTime = System.nanoTime();
    Complex.dftBlocked(dataBlocked, dft);
    dftTime = System.nanoTime() - dftTime;
    String dftString = " DFT Time: " + dftTime * 1.0e-9 + " s\n";

    long fftTime = System.nanoTime();
    double[] data2 = new double[2 * n * nFFTs];
    int index = 0;
    for (int ii = 0; ii < n; ii += stride) {
      for (int i = 0; i < nFFTs; i++) {
        data2[index++] = dataBlocked[ii];
      }
    }

    complex.fft(data2, offset, stride, nextFFT);
    fftTime = System.nanoTime() - fftTime;
    String fftString = " FFT Time: " + fftTime * 1.0e-9 + " s";

    // Test the FFT is equals the DFT result.
    index = 0;
    for (int ii = 0; ii < n; ii++) {
      for (int i = 0; i < nFFTs; i++) {
        assertEquals(" FFT " + i + " Forward " + info + " at position: " + ii,
            dft[ii], data2[index], tolerance);
        assertEquals(" FFT " + i + " Forward " + info + " at position: " + ii,
            dft[ii + n], data2[index + im], tolerance);
        index++;
      }
    }

    // The FFT is faster than the DFT.
    String message = fftString + dftString;
    // assertTrue(message, fftTime < dftTime);

    // Test that X = IFFT(FFT(X)).
    complex.inverse(data2, offset, stride, nextFFT);

    // Test the FFT is equals the DFT result.
    index = 0;
    for (int ii = 0; ii < nextFFT; ii += 2) {
      for (int i = 0; i < nFFTs; i++) {
        assertEquals(" IFFT(FFT(X)) " + info + " at position: " + ii + " FFT " + i,
            orig[ii], data2[index], tolerance);
        index++;
      }
    }
  }

}
