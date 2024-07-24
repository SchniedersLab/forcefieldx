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

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import javax.annotation.Nullable;
import java.util.Arrays;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.util.Objects.requireNonNullElseGet;

/**
 * Compute the 3D FFT of complex, double precision input of arbitrary dimensions via 1D Mixed Radix
 * FFTs in parallel.
 *
 * <p>The location of the input point [i, j, k] within the input array must be: <br>
 * double real = input[x*nextX + y*nextY + z*nextZ] <br>
 * double imag = input[x*nextX + y*nextY + z*nextZ + 1] <br>
 * where <br>
 * int nextX = 2 <br>
 * int nextY = 2*nX <br>
 * int nextZ = 2*nX*nY <br>
 *
 * @author Michal J. Schnieders
 * @see Complex
 * @since 1.0
 */
public class Complex3DParallel {

  private static final Logger logger = Logger.getLogger(Complex3DParallel.class.getName());
  /**
   * The X-dimension.
   */
  private final int nX;
  /**
   * The Y-dimension.
   */
  private final int nY;
  /**
   * The Z-dimension.
   */
  private final int nZ;
  /**
   * The offset from any real value to its imaginary part.
   */
  private final int im;
  /**
   * The offset from any real value to the next real value (1 for blocked and 2 for interleaved).
   */
  private final int ii;
  /**
   * The offset to the next X value.
   */
  private final int nextX;
  /**
   * The offset to the next Y value.
   */
  private final int nextY;
  /**
   * The offset to the next Z value.
   */
  private final int nextZ;
  /**
   * The offset to the next X value for the transposed data.
   */
  private final int trNextX;
  /**
   * The offset to the next Y value for the transposed data.
   */
  private final int trNextY;
  /**
   * The offset to the next Z value for the transposed data.
   */
  private final int trNextZ;
  /**
   * The reciprocal space data.
   */
  private final double[] recip;
  /**
   * The number of threads to use.
   */
  private final int threadCount;
  /**
   * The ParallelTeam instance.
   */
  private final ParallelTeam parallelTeam;
  /**
   * The 2D FFTs for each thread for XY-planes.
   */
  private final Complex2D[] fftXY;
  /**
   * The 1D FFTs for each thread to compute FFTs along the Z-dimension.
   * Each thread operates one YZ-plane at a time to complete nY FFTs.
   */
  private final Complex[] fftZ;
  /**
   * The FFTs along the Z-dimension will use offset between real and imaginary parts.
   */
  private final int internalImZ;
  /**
   * The IntegerSchedule to use.
   */
  private final IntegerSchedule schedule;
  /**
   * The X-dimension minus 1.
   */
  private final int nXm1;
  /**
   * The Y-dimension minus 1.
   */
  private final int nYm1;
  /**
   * The Z-dimension minus 1.
   */
  private final int nZm1;
  /**
   * The FFTRegion for parallel calculation of a 3D FFT.
   */
  private final FFTRegion fftRegion;
  /**
   * The IFFTRegion for parallel calculation of an inverse 3D FFT.
   */
  private final IFFTRegion ifftRegion;
  /**
   * The ConvolutionRegion for parallel calculation of a convolution.
   */
  private final ConvolutionRegion convRegion;
  /**
   * This is a reference to the input array for convenience.
   * The input array must be of size 2 * nX * nY * nZ.
   */
  public double[] input;
  /**
   * Use SIMD instructions if available.
   */
  private boolean useSIMD;
  /**
   * Pack the FFTs for optimal use of SIMD instructions.
   */
  private boolean packFFTs;
  /**
   * Work array used for transforms along the Z-axis.
   */
  public final double[] work3D;

  /**
   * Initialize the 3D FFT for complex 3D matrix.
   *
   * @param nX           X-dimension.
   * @param nY           Y-dimension.
   * @param nZ           Z-dimension.
   * @param parallelTeam A ParallelTeam instance.
   * @since 1.0
   */
  public Complex3DParallel(int nX, int nY, int nZ, ParallelTeam parallelTeam) {
    this(nX, nY, nZ, parallelTeam, DataLayout3D.INTERLEAVED);
  }

  /**
   * Initialize the 3D FFT for complex 3D matrix.
   *
   * @param nX           X-dimension.
   * @param nY           Y-dimension.
   * @param nZ           Z-dimension.
   * @param parallelTeam A ParallelTeam instance.
   * @param dataLayout   The data layout.
   * @since 1.0
   */
  public Complex3DParallel(int nX, int nY, int nZ, ParallelTeam parallelTeam, DataLayout3D dataLayout) {
    this(nX, nY, nZ, parallelTeam, null, dataLayout);
  }

  /**
   * Initialize the 3D FFT for complex 3D matrix.
   *
   * @param nX              X-dimension.
   * @param nY              Y-dimension.
   * @param nZ              Z-dimension.
   * @param parallelTeam    A ParallelTeam instance.
   * @param integerSchedule The IntegerSchedule to use.
   * @since 1.0
   */
  public Complex3DParallel(int nX, int nY, int nZ, ParallelTeam parallelTeam, @Nullable IntegerSchedule integerSchedule) {
    this(nX, nY, nZ, parallelTeam, integerSchedule, DataLayout3D.INTERLEAVED);
  }

  /**
   * Initialize the 3D FFT for complex 3D matrix.
   *
   * @param nX              X-dimension.
   * @param nY              Y-dimension.
   * @param nZ              Z-dimension.
   * @param parallelTeam    A ParallelTeam instance.
   * @param integerSchedule The IntegerSchedule to use.
   * @param dataLayout      The data layout.
   * @since 1.0
   */
  public Complex3DParallel(int nX, int nY, int nZ, ParallelTeam parallelTeam,
                           @Nullable IntegerSchedule integerSchedule, DataLayout3D dataLayout) {
    this.nX = nX;
    this.nY = nY;
    this.nZ = nZ;
    this.parallelTeam = parallelTeam;
    recip = new double[nX * nY * nZ];

    DataLayout1D dataLayout1D;
    DataLayout2D dataLayout2D;
    switch (dataLayout) {
      default:
      case INTERLEAVED:
        // Interleaved data layout.
        im = 1;
        ii = 2;
        nextX = 2;
        nextY = 2 * nX;
        nextZ = 2 * nX * nY;
        // Internal FFTs will be performed in interleaved format.
        dataLayout1D = DataLayout1D.INTERLEAVED;
        dataLayout2D = DataLayout2D.INTERLEAVED;
        // Transforms along the Z-axis will be repacked into 1D interleaved format.
        internalImZ = 1;
        trNextY = 2;
        trNextZ = 2 * nY;
        trNextX = 2 * nY * nZ;
        break;
      case BLOCKED_X:
        // Blocking is along the X-axis.
        im = nX;
        ii = 1;
        nextX = 1;
        nextY = 2 * nX;
        nextZ = 2 * nX * nY;
        // Internal FFTs will be performed in blocked format.
        dataLayout1D = DataLayout1D.BLOCKED;
        dataLayout2D = DataLayout2D.BLOCKED_X;
        // Transforms along the Z-axis will be repacked into 1D blocked format.
        internalImZ = nZ * nY;
        trNextY = 1;
        trNextZ = nY;
        trNextX = 2 * nY * nZ;
        break;
      case BLOCKED_XY:
        // Blocking is based on 2D XY-planes.
        im = nX * nY;
        ii = 1;
        nextX = 1;
        nextY = nX;
        nextZ = 2 * nY * nX;
        // Internal FFTs will be performed in blocked format.
        dataLayout1D = DataLayout1D.BLOCKED;
        dataLayout2D = DataLayout2D.BLOCKED_XY;
        // Transforms along the Z-axis will be repacked into 1D blocked format.
        internalImZ = nY * nZ;
        trNextY = 1;
        trNextZ = nY;
        trNextX = 2 * nY * nZ;
        break;
      case BLOCKED_XYZ:
        // Blocking is based on 3D XYZ-volume with all real values followed by all imaginary.
        im = nX * nY * nZ;
        ii = 1;
        nextX = 1;
        nextY = nX;
        nextZ = nY * nX;
        // Internal 1D FFTs will be performed in blocked format.
        dataLayout1D = DataLayout1D.BLOCKED;
        dataLayout2D = DataLayout2D.BLOCKED_XY;
        // Transforms along the Z-axis will be repacked into 1D blocked format.
        internalImZ = nY * nZ;
        trNextY = 1;
        trNextZ = nY;
        trNextX = 2 * nY * nZ;
        break;
    }

    nXm1 = this.nX - 1;
    nYm1 = this.nY - 1;
    nZm1 = this.nZ - 1;
    threadCount = parallelTeam.getThreadCount();
    schedule = requireNonNullElseGet(integerSchedule, IntegerSchedule::fixed);

    // Do not use SIMD by default for now.
    useSIMD = false;
    String simd = System.getProperty("fft.useSIMD", Boolean.toString(useSIMD));
    try {
      useSIMD = Boolean.parseBoolean(simd);
    } catch (Exception e) {
      useSIMD = false;
    }

    packFFTs = false;
    String pack = System.getProperty("fft.packFFTs", Boolean.toString(packFFTs));
    try {
      packFFTs = Boolean.parseBoolean(pack);
    } catch (Exception e) {
      packFFTs = false;
    }

    fftXY = new Complex2D[threadCount];
    for (int i = 0; i < threadCount; i++) {
      fftXY[i] = new Complex2D(nX, nY, dataLayout2D, im);
      fftXY[i].setPackFFTs(packFFTs);
      fftXY[i].setUseSIMD(useSIMD);
    }

    fftZ = new Complex[threadCount];
    for (int i = 0; i < threadCount; i++) {
      fftZ[i] = new Complex(nZ, dataLayout1D, internalImZ, nY);
      fftZ[i].setUseSIMD(useSIMD);
    }
    work3D = new double[2 * nX * nY * nZ];
    fftRegion = new FFTRegion();
    ifftRegion = new IFFTRegion();
    convRegion = new ConvolutionRegion();
  }

  public String toString() {
    return "Complex3DParallel {" +
        "nX=" + nX +
        ", nY=" + nY +
        ", nZ=" + nZ +
        ", im=" + im +
        ", ii=" + ii +
        ", nextX=" + nextX +
        ", nextY=" + nextY +
        ", nextZ=" + nextZ +
        ", threadCount=" + threadCount +
        ", parallelTeam=" + parallelTeam +
        ", internalImZ=" + internalImZ +
        ", schedule=" + schedule +
        ", nXm1=" + nXm1 +
        ", nYm1=" + nYm1 +
        ", nZm1=" + nZm1 +
        ", fftRegion=" + fftRegion +
        ", ifftRegion=" + ifftRegion +
        ", convRegion=" + convRegion +
        ", useSIMD=" + useSIMD +
        ", packFFTs=" + packFFTs +
        '}';
  }

  /**
   * Use SIMD instructions if available.
   *
   * @param useSIMD True to use SIMD instructions.
   */
  public void setUseSIMD(boolean useSIMD) {
    this.useSIMD = useSIMD;
    for (int i = 0; i < threadCount; i++) {
      fftXY[i].setUseSIMD(useSIMD);
      fftZ[i].setUseSIMD(useSIMD);
    }
  }

  /**
   * Pack the FFTs for optimal use of SIMD instructions.
   *
   * @param packFFTs True to pack the FFTs.
   */
  public void setPackFFTs(boolean packFFTs) {
    this.packFFTs = packFFTs;
    for (int i = 0; i < threadCount; i++) {
      fftXY[i].setPackFFTs(packFFTs);
    }
  }

  /**
   * Compute the 3D FFT, perform a multiplication in reciprocal space,
   * and the inverse 3D FFT in parallel.
   *
   * @param input The input array must be of size 2 * nX * nY * nZ.
   * @since 1.0
   */
  public void convolution(final double[] input) {
    this.input = input;
    try {
      parallelTeam.execute(convRegion);
    } catch (Exception e) {
      String message = "Fatal exception evaluating a convolution.\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  /**
   * Compute the 3D FFT in parallel.
   *
   * @param input The input array must be of size 2 * nX * nY * nZ.
   * @since 1.0
   */
  public void fft(final double[] input) {
    this.input = input;
    try {
      parallelTeam.execute(fftRegion);
    } catch (Exception e) {
      String message = " Fatal exception evaluating the FFT.\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  /**
   * Get the timings for each thread.
   *
   * @return The timings for each thread.
   */
  public long[] getTiming() {
    return convRegion.getTiming();
  }

  /**
   * Compute the inverse 3D FFT in parallel.
   *
   * @param input The input array must be of size 2 * nX * nY * nZ.
   * @since 1.0
   */
  public void ifft(final double[] input) {
    this.input = input;
    try {
      parallelTeam.execute(ifftRegion);
    } catch (Exception e) {
      String message = "Fatal exception evaluating the inverse FFT.\n";
      logger.log(Level.SEVERE, message, e);
      System.exit(-1);
    }
  }

  /**
   * Initialize the timing array.
   */
  public void initTiming() {
    convRegion.initTiming();
  }

  public String timingString() {
    return convRegion.timingString();
  }

  /**
   * Setter for the field <code>recip</code>.
   *
   * @param recip an array of double.
   */
  public void setRecip(double[] recip) {
    // Reorder the reciprocal space data for convolution.
    // Input
    // trNextY = ii
    // trNextZ = nY*ii
    // real[y, z] = work[y*trNextY + z*trNextZ]
    // imag[y, z] = work[y*trNextY + z*trNextZ + internalImZ]
    int recipNextY = nX;
    int recipNextZ = nY * nX;
    int index = 0;
    for (int x = 0; x < nX; x++) {
      int dx = x;
      for (int z = 0; z < nZ; z++) {
        int dz = dx + z * recipNextZ;
        for (int y = 0; y < nY; y++) {
          int conv = y * recipNextY + dz;
          this.recip[index] = recip[conv];
          index++;
        }
      }
    }
  }

  /**
   * An external ParallelRegion can be used as follows: <code>
   * start() {
   * fftRegion.input = input;
   * }
   * run(){
   * execute(0, nZm1, fftRegion.fftXYLoop[threadID]);
   * execute(0, nXm1, fftRegion.fftZLoop[threadID]);
   * }
   * </code>
   */
  private class FFTRegion extends ParallelRegion {

    private final FFTXYLoop[] fftXYLoop;
    private final FFTZLoop[] fftZLoop;
    private final TransposeLoop[] transposeLoop;
    private final UnTransposeLoop[] unTransposeLoop;

    private FFTRegion() {
      fftXYLoop = new FFTXYLoop[threadCount];
      fftZLoop = new FFTZLoop[threadCount];
      transposeLoop = new TransposeLoop[threadCount];
      unTransposeLoop = new UnTransposeLoop[threadCount];
      for (int i = 0; i < threadCount; i++) {
        fftXYLoop[i] = new FFTXYLoop();
        fftZLoop[i] = new FFTZLoop();
        transposeLoop[i] = new TransposeLoop();
        unTransposeLoop[i] = new UnTransposeLoop();
      }
    }

    @Override
    public void run() {
      int threadIndex = getThreadIndex();
      try {
        execute(0, nZm1, fftXYLoop[threadIndex]);
        execute(0, nXm1, transposeLoop[threadIndex]);
        execute(0, nXm1, fftZLoop[threadIndex]);
        execute(0, nZm1, unTransposeLoop[threadIndex]);
      } catch (Exception e) {
        logger.severe(e.toString());
      }
    }
  }

  /**
   * An external ParallelRegion can be used as follows: <code>
   * start() {
   * ifftRegion.input = input;
   * }
   * run(){
   * execute(0, nXm1, ifftRegion.ifftZLoop[threadID]);
   * execute(0, nZm1, ifftRegion.ifftXYLoop[threadID]);
   * }
   * </code>
   */
  private class IFFTRegion extends ParallelRegion {

    private final IFFTXYLoop[] ifftXYLoop;
    private final IFFTZLoop[] ifftZLoop;
    private final TransposeLoop[] transposeLoop;
    private final UnTransposeLoop[] unTransposeLoop;

    private IFFTRegion() {
      ifftXYLoop = new IFFTXYLoop[threadCount];
      transposeLoop = new TransposeLoop[threadCount];
      ifftZLoop = new IFFTZLoop[threadCount];
      unTransposeLoop = new UnTransposeLoop[threadCount];
      for (int i = 0; i < threadCount; i++) {
        ifftXYLoop[i] = new IFFTXYLoop();
        transposeLoop[i] = new TransposeLoop();
        ifftZLoop[i] = new IFFTZLoop();
        unTransposeLoop[i] = new UnTransposeLoop();
      }
    }

    @Override
    public void run() {
      int threadIndex = getThreadIndex();
      try {
        execute(0, nXm1, transposeLoop[threadIndex]);
        execute(0, nXm1, ifftZLoop[threadIndex]);
        execute(0, nZm1, unTransposeLoop[threadIndex]);
        execute(0, nZm1, ifftXYLoop[threadIndex]);
      } catch (Exception e) {
        logger.severe(e.toString());
      }
    }
  }

  /**
   * An external ParallelRegion can be used as follows: <code>
   * start() {
   * convRegion.input = input;
   * }
   * run(){
   * execute(0, nZm1, convRegion.fftXYLoop[threadID]);
   * execute(0, nYm1, convRegion.fftZIZLoop[threadID]);
   * execute(0, nZm1, convRegion.ifftXYLoop[threadID]);
   * }
   * </code>
   */
  private class ConvolutionRegion extends ParallelRegion {

    private final FFTXYLoop[] fftXYLoop;
    private final TransposeLoop[] transposeLoop;
    private final FFTZIZLoop[] fftZIZLoop;
    private final UnTransposeLoop[] unTransposeLoop;
    private final IFFTXYLoop[] ifftXYLoop;
    private final long[] convTime;

    private ConvolutionRegion() {
      fftXYLoop = new FFTXYLoop[threadCount];
      transposeLoop = new TransposeLoop[threadCount];
      fftZIZLoop = new FFTZIZLoop[threadCount];
      unTransposeLoop = new UnTransposeLoop[threadCount];
      ifftXYLoop = new IFFTXYLoop[threadCount];
      convTime = new long[threadCount];
      for (int i = 0; i < threadCount; i++) {
        fftXYLoop[i] = new FFTXYLoop();
        transposeLoop[i] = new TransposeLoop();
        fftZIZLoop[i] = new FFTZIZLoop();
        unTransposeLoop[i] = new UnTransposeLoop();
        ifftXYLoop[i] = new IFFTXYLoop();

      }
    }

    public void initTiming() {
      for (int i = 0; i < threadCount; i++) {
        fftXYLoop[i].time = 0;
        transposeLoop[i].time = 0;
        fftZIZLoop[i].time = 0;
        unTransposeLoop[i].time = 0;
        ifftXYLoop[i].time = 0;
      }
    }

    public long[] getTiming() {
      for (int i = 0; i < threadCount; i++) {
        convTime[i] = convRegion.fftXYLoop[i].time
            + convRegion.transposeLoop[i].time
            + convRegion.fftZIZLoop[i].time
            + convRegion.unTransposeLoop[i].time
            + convRegion.ifftXYLoop[i].time;
      }
      return convTime;
    }

    public String timingString() {
      StringBuilder sb = new StringBuilder();
      double xysum = 0.0;
      double transsum = 0.0;
      double zizsum = 0.0;
      double untranssum = 0.0;
      double ixysum = 0.0;
      for (int i = 0; i < threadCount; i++) {
        double fftxy = fftXYLoop[i].getTime() * 1e-9;
        double trans = transposeLoop[i].getTime() * 1e-9;
        double ziz = fftZIZLoop[i].getTime() * 1e-9;
        double untrans = unTransposeLoop[i].getTime() * 1e-9;
        double ifftxy = ifftXYLoop[i].getTime() * 1e-9;
        String s = String.format(
            "  Thread %3d: FFTXY=%8.6f, Trans=%8.6f, FFTZIZ=%8.6f, UnTrans=%8.6f, IFFTXY=%8.6f\n",
            i, fftxy, trans, ziz, untrans, ifftxy);
        sb.append(s);
        xysum += fftxy;
        transsum += trans;
        zizsum += ziz;
        untranssum += untrans;
        ixysum += ifftxy;
      }
      String s = String.format(
          "  Sum       : FFTXY=%8.6f, Trans=%8.6f, FFTZIZ=%8.6f, UnTrans=%8.6f, IFFTXY=%8.6f\n",
          xysum, transsum, zizsum, untranssum, ixysum);
      sb.append(s);
      return sb.toString();
    }

    @Override
    public void run() {
      int threadIndex = getThreadIndex();
      try {
        execute(0, nZm1, fftXYLoop[threadIndex]);
        execute(0, nXm1, transposeLoop[threadIndex]);
        execute(0, nXm1, fftZIZLoop[threadIndex]);
        execute(0, nZm1, unTransposeLoop[threadIndex]);
        execute(0, nZm1, ifftXYLoop[threadIndex]);
      } catch (Exception e) {
        logger.severe(e.toString());
      }
    }
  }

  private class FFTXYLoop extends IntegerForLoop {

    private Complex2D localFFTXY;
    private long time;

    @Override
    public void run(final int lb, final int ub) {
      for (int z = lb; z <= ub; z++) {
        int offset = z * nextZ;
        localFFTXY.fft(input, offset);
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return schedule;
    }

    public long getTime() {
      return time;
    }

    public void finish() {
      time += System.nanoTime();
    }

    @Override
    public void start() {
      time -= System.nanoTime();
      localFFTXY = fftXY[getThreadIndex()];
    }
  }

  private class IFFTXYLoop extends IntegerForLoop {

    private Complex2D localFFTXY;
    private long time;

    @Override
    public void run(final int lb, final int ub) {
      for (int z = lb; z <= ub; z++) {
        int offset = z * nextZ;
        localFFTXY.ifft(input, offset);
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return schedule;
    }

    public long getTime() {
      return time;
    }

    public void finish() {
      time += System.nanoTime();
    }

    @Override
    public void start() {
      time -= System.nanoTime();
      localFFTXY = fftXY[getThreadIndex()];
    }
  }

  private class FFTZLoop extends IntegerForLoop {

    private Complex localFFTZ;
    private long time;

    private FFTZLoop() {
      // Empty.
    }

    @Override
    public void run(final int lb, final int ub) {
      for (int x = lb; x <= ub; x++) {
        int offset = x * nY * nZ * ii;
        localFFTZ.fft(work3D, offset, ii);
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return schedule;
    }

    public long getTime() {
      return time;
    }

    public void finish() {
      time += System.nanoTime();
    }

    @Override
    public void start() {
      time -= System.nanoTime();
      int threadID = getThreadIndex();
      localFFTZ = fftZ[threadID];
    }
  }

  private class IFFTZLoop extends IntegerForLoop {

    private Complex localFFTZ;

    private IFFTZLoop() {
      // Empty.
    }

    @Override
    public void run(final int lb, final int ub) {
      for (int x = lb; x <= ub; x++) {
        int offset = x * nY * nZ * ii;
        localFFTZ.ifft(work3D, offset, ii);
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return schedule;
    }

    @Override
    public void start() {
      int threadID = getThreadIndex();
      localFFTZ = fftZ[threadID];
    }
  }

  private class FFTZIZLoop extends IntegerForLoop {

    private Complex localFFTZ;
    private long time;

    private FFTZIZLoop() {
      // Empty.
    }

    @Override
    public void run(final int lb, final int ub) {
      for (int x = lb; x <= ub; x++) {
        // Offset into the work3D array.
        int offset = x * nY * nZ * ii;
        localFFTZ.fft(work3D, offset, ii);
        recipConv(x, work3D);
        localFFTZ.ifft(work3D, offset, ii);
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return schedule;
    }

    public long getTime() {
      return time;
    }

    public void finish() {
      time += System.nanoTime();
    }

    @Override
    public void start() {
      time -= System.nanoTime();
      int threadID = getThreadIndex();
      localFFTZ = fftZ[threadID];
    }
  }

  private class TransposeLoop extends IntegerForLoop {

    private long time;

    @Override
    public void run(final int lb, final int ub) {
      // The input data order is X,Y,Z (X varies most quickly)
      // real[x, y, z] = input[x*nextX + y*nextY + z*nextZ]
      // imag[x, y, z] = input[x*nextX + y*nextY + z*nextZ + im]
      // The output data order is Y,Z,X (Y varies most quickly)
      // real[x, y, z] = input[y*trNextY + z*trNextZ + x*trNextX]
      // imag[x, y, z] = input[y*trNextY + z*trNextZ + x*trNextX + im]

      // Starting output index.
      int index = lb * nY * nZ * ii;
      for (int x = lb; x <= ub; x++) {
        int dx = x * nextX;
        for (int z = 0; z < nZ; z++) {
          int dz = dx + z * nextZ;
          for (int y = 0; y < nY; y++) {
            double real = input[y * nextY + dz];
            double imag = input[y * nextY + dz + im];
            work3D[index] = real;
            work3D[index + internalImZ] = imag;
            index += ii;
          }
        }
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return schedule;
    }

    public long getTime() {
      return time;
    }

    @Override
    public void finish() {
      time += System.nanoTime();
    }

    @Override
    public void start() {
      time -= System.nanoTime();
    }
  }

  private class UnTransposeLoop extends IntegerForLoop {

    private long time;

    @Override
    public void run(final int lb, final int ub) {
      // The transposed input data order is Y,Z,X (Y varies most quickly)
      // real[x, y, z] = input[y*trNextY + z*trNextZ + x*trNextX]
      // imag[x, y, z] = input[y*trNextY + z*trNextZ + x*trNextX + im]
      // The output data order is X,Y,Z (X varies most quickly)
      // real[x, y, z] = input[x*nextX + y*nextY + z*nextZ]
      // imag[x, y, z] = input[x*nextX + y*nextY + z*nextZ + im]
      int index = lb * nextZ;
      for (int z = lb; z <= ub; z++) {
        int dz = z * trNextZ;
        for (int y = 0; y < nY; y++) {
          int dy = y * trNextY + dz;
          for (int x = 0; x < nX; x++) {
            int dx = x * trNextX + dy;
            double real = work3D[dx];
            double imag = work3D[dx + internalImZ];
            input[index] = real;
            input[index + im] = imag;
            index += ii;
          }
        }
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return schedule;
    }

    public long getTime() {
      return time;
    }

    @Override
    public void finish() {
      time += System.nanoTime();
    }

    @Override
    public void start() {
      time -= System.nanoTime();
    }
  }
  
  /**
   * Perform a multiplication by the reciprocal space data.
   *
   * @param x    The X-value for this Y-Z plane.
   * @param work The input array.
   */
  private void recipConv(int x, double[] work) {
    int index = x * nY * nZ * ii;
    int rindex = x * nY * nZ;
    for (int i = 0; i < nY * nZ; i++) {
      double r = recip[rindex++];
      work[index] *= r;
      work[index + internalImZ] *= r;
      index += ii;
    }
  }

  /**
   * Initialize a 3D data for testing purposes.
   *
   * @param dim The dimension of the cube.
   * @since 1.0
   */
  public static double[] initRandomData(int dim, ParallelTeam parallelTeam) {
    int n = dim * dim * dim;
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
                        int index = dim * dim * lb * 2;
                        for (int i = lb; i <= ub; i++) {
                          for (int j = 0; j < dim; j++) {
                            for (int k = 0; k < dim; k++) {
                              double randomNumber = randomNumberGenerator.nextDouble();
                              data[index] = randomNumber;
                              index += 2;
                            }
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

  /**
   * Test the Complex3DParallel FFT.
   *
   * @param args an array of {@link java.lang.String} objects.
   * @throws java.lang.Exception if any.
   * @since 1.0
   */
  public static void main(String[] args) throws Exception {
    int dimNotFinal = 128;
    int nCPU = ParallelTeam.getDefaultThreadCount();
    int reps = 5;
    boolean blocked = false;
    try {
      dimNotFinal = Integer.parseInt(args[0]);
      if (dimNotFinal < 1) {
        dimNotFinal = 100;
      }
      nCPU = Integer.parseInt(args[1]);
      if (nCPU < 1) {
        nCPU = ParallelTeam.getDefaultThreadCount();
      }
      reps = Integer.parseInt(args[2]);
      if (reps < 1) {
        reps = 5;
      }
      blocked = Boolean.parseBoolean(args[3]);
    } catch (Exception e) {
      //
    }
    final int dim = dimNotFinal;
    System.out.printf("Initializing a %d cubed grid for %d CPUs.\n"
            + "The best timing out of %d repetitions will be used.%n",
        dim, nCPU, reps);
    // One dimension of the serial array divided by the number of threads.
    Complex3D complex3D;
    Complex3DParallel complex3DParallel;
    ParallelTeam parallelTeam = new ParallelTeam(nCPU);
    if (blocked) {
      complex3D = new Complex3D(dim, dim, dim, DataLayout3D.BLOCKED_X);
      complex3DParallel = new Complex3DParallel(dim, dim, dim, parallelTeam, DataLayout3D.BLOCKED_X);
    } else {
      complex3D = new Complex3D(dim, dim, dim, DataLayout3D.INTERLEAVED);
      complex3DParallel = new Complex3DParallel(dim, dim, dim, parallelTeam, DataLayout3D.INTERLEAVED);
    }
    final int dimCubed = dim * dim * dim;
    final double[] data = initRandomData(dim, parallelTeam);
    final double[] work = new double[dimCubed];
    Arrays.fill(work, 1.0);

    double toSeconds = 0.000000001;
    long seqTime = Long.MAX_VALUE;
    long parTime = Long.MAX_VALUE;
    long seqTimeConv = Long.MAX_VALUE;
    long parTimeConv = Long.MAX_VALUE;

    complex3D.setRecip(work);
    complex3DParallel.setRecip(work);

    // Warm-up
    System.out.println("Warm Up Sequential FFT");
    complex3D.fft(data);
    System.out.println("Warm Up Sequential IFFT");
    complex3D.ifft(data);
    System.out.println("Warm Up Sequential Convolution");
    complex3D.convolution(data);
    for (int i = 0; i < reps; i++) {
      System.out.printf(" Iteration %d%n", i + 1);
      long time = System.nanoTime();
      complex3D.fft(data);
      complex3D.ifft(data);
      time = (System.nanoTime() - time);
      System.out.printf("  Sequential FFT:  %9.6f (sec)%n", toSeconds * time);
      if (time < seqTime) {
        seqTime = time;
      }
      time = System.nanoTime();
      complex3D.convolution(data);
      time = (System.nanoTime() - time);
      System.out.printf("  Sequential Conv: %9.6f (sec)%n", toSeconds * time);
      if (time < seqTimeConv) {
        seqTimeConv = time;
      }
    }

    // Warm-up
    System.out.println("Warm up Parallel FFT");
    complex3DParallel.fft(data);
    System.out.println("Warm up Parallel IFFT");
    complex3DParallel.ifft(data);
    System.out.println("Warm up Parallel Convolution");
    complex3DParallel.convolution(data);
    complex3DParallel.initTiming();

    for (int i = 0; i < reps; i++) {
      System.out.printf(" Iteration %d%n", i + 1);
      long time = System.nanoTime();
      complex3DParallel.fft(data);
      complex3DParallel.ifft(data);
      time = (System.nanoTime() - time);
      System.out.printf("  Parallel FFT:  %9.6f (sec)%n", toSeconds * time);
      if (time < parTime) {
        parTime = time;
      }

      time = System.nanoTime();
      complex3DParallel.convolution(data);
      time = (System.nanoTime() - time);
      System.out.printf("  Parallel Conv: %9.6f (sec)%n", toSeconds * time);
      if (time < parTimeConv) {
        parTimeConv = time;
      }
    }

    System.out.printf(" Best Sequential FFT Time:   %9.6f (sec)%n", toSeconds * seqTime);
    System.out.printf(" Best Sequential Conv. Time: %9.6f (sec)%n", toSeconds * seqTimeConv);
    System.out.printf(" Best Parallel FFT Time:     %9.6f (sec)%n", toSeconds * parTime);
    System.out.printf(" Best Parallel Conv. Time:   %9.6f (sec)%n", toSeconds * parTimeConv);
    System.out.printf(" 3D FFT Speedup:             %9.6f X%n", (double) seqTime / parTime);
    System.out.printf(" 3D Conv Speedup:            %9.6f X%n", (double) seqTimeConv / parTimeConv);

    System.out.printf(" Parallel Convolution Timings:\n" + complex3DParallel.timingString());

    parallelTeam.shutdown();
  }
}
