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
 * Compare the forward Complex1D FFT against the Complex FFT.
 *
 * @author Michael J. Schnieders
 */
@RunWith(Parameterized.class)
public class Complex1DTest extends FFXTest {

  private final String info;
  private final int n;

  public Complex1DTest(String info, int n) {
    this.info = info;
    this.n = n;
  }

  @Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][]{
            {"Test n = 6", 6},
            {"Test n = 8", 8},
            {"Test n = 10", 10},
            {"Test n = 12", 12},
            {"Test n = 15", 15},
            {"Test n = 18", 18},
            {"Test n = 20", 20},
            {"Test n = 21", 21},
            {"Test n = 27", 27},
            {"Test n = 30", 30},
            {"Test n = 32", 32},
            {"Test n = 38", 38},
            {"Test n = 45", 45}
        });
  }

  @Test
  public void testForwardFFTMatchesComplex() {
    double tolerance = 1.0e-11;
    double[] expected = createRandomComplexData();
    double[] actual = expected.clone();

    Complex complex = new Complex(n);
    Complex1D complex1D = new Complex1D(n);
    complex.fft(expected, 0, 2);
    complex1D.fft(actual, 0, 2);

    for (int i = 0; i < 2 * n; i++) {
      assertEquals(info + " at position: " + i, expected[i], actual[i], tolerance);
    }
  }

  @Test
  public void testInverseFFTMatchesComplex() {
    double tolerance = 1.0e-11;
    double[] expected = createRandomComplexData();
    double[] actual = expected.clone();

    Complex complex = new Complex(n);
    Complex1D complex1D = new Complex1D(n);
    complex.ifft(expected, 0, 2);
    complex1D.ifft(actual, 0, 2);

    for (int i = 0; i < 2 * n; i++) {
      assertEquals(info + " inverse at position: " + i, expected[i], actual[i], tolerance);
    }
  }

  @Test
  public void testForwardInverseRoundTrip() {
    double tolerance = 1.0e-11;
    double[] original = createRandomComplexData();
    double[] actual = original.clone();

    Complex1D complex1D = new Complex1D(n);
    complex1D.fft(actual, 0, 2);
    complex1D.ifft(actual, 0, 2);

    for (int i = 0; i < 2 * n; i++) {
      assertEquals(info + " round-trip at position: " + i, original[i], actual[i] / n, tolerance);
    }
  }

  @Test
  public void testBlockedForwardFFTMatchesComplex() {
    if (!Complex.preferredDimension(n)) {
      return;
    }
    double tolerance = 1.0e-11;
    double[] expected = createRandomBlockedComplexData();
    double[] actual = expected.clone();

    Complex complex = new Complex(n, DataLayout1D.BLOCKED, n);
    Complex1D complex1D = new Complex1D(n, DataLayout1D.BLOCKED, n);
    complex.fft(expected, 0, 1);
    complex1D.fft(actual, 0, 1);

    for (int i = 0; i < 2 * n; i++) {
      assertEquals(info + " blocked at position: " + i, expected[i], actual[i], tolerance);
    }
  }

  @Test
  public void testBlockedInverseFFTMatchesComplex() {
    if (!Complex.preferredDimension(n)) {
      return;
    }
    double tolerance = 1.0e-11;
    double[] expected = createRandomBlockedComplexData();
    double[] actual = expected.clone();

    Complex complex = new Complex(n, DataLayout1D.BLOCKED, n);
    Complex1D complex1D = new Complex1D(n, DataLayout1D.BLOCKED, n);
    complex.ifft(expected, 0, 1);
    complex1D.ifft(actual, 0, 1);

    for (int i = 0; i < 2 * n; i++) {
      assertEquals(info + " blocked inverse at position: " + i, expected[i], actual[i], tolerance);
    }
  }

  @Test
  public void testBlockedForwardInverseRoundTrip() {
    if (!Complex.preferredDimension(n)) {
      return;
    }
    double tolerance = 1.0e-11;
    double[] original = createRandomBlockedComplexData();
    double[] actual = original.clone();

    Complex1D complex1D = new Complex1D(n, DataLayout1D.BLOCKED, n);
    complex1D.fft(actual, 0, 1);
    complex1D.ifft(actual, 0, 1);

    for (int i = 0; i < 2 * n; i++) {
      assertEquals(info + " blocked round-trip at position: " + i, original[i], actual[i] / n, tolerance);
    }
  }

  private double[] createRandomComplexData() {
    double[] data = new double[2 * n];
    Random random = new Random(1);
    for (int i = 0; i < n; i++) {
      int index = 2 * i;
      data[index] = random.nextDouble();
      data[index + 1] = random.nextDouble();
    }
    return data;
  }

  private double[] createRandomBlockedComplexData() {
    double[] data = new double[2 * n];
    Random random = new Random(1);
    for (int i = 0; i < n; i++) {
      data[i] = random.nextDouble();
      data[i + n] = random.nextDouble();
    }
    return data;
  }
}
