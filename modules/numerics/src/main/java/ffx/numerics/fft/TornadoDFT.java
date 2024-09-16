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

import ffx.numerics.tornado.FFXTornado;
import uk.ac.manchester.tornado.api.ImmutableTaskGraph;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.annotations.Parallel;
import uk.ac.manchester.tornado.api.common.TornadoDevice;
import uk.ac.manchester.tornado.api.runtime.TornadoRuntimeProvider;

import static uk.ac.manchester.tornado.api.enums.DataTransferMode.EVERY_EXECUTION;
import static uk.ac.manchester.tornado.api.math.TornadoMath.cos;
import static uk.ac.manchester.tornado.api.math.TornadoMath.floatPI;
import static uk.ac.manchester.tornado.api.math.TornadoMath.sin;

/**
 * Proof-of-concept use of the TornadoVM for parallelization of Java code.
 */
public class TornadoDFT {

  float[] inReal;
  float[] inImag;
  float[] outReal;
  float[] outImag;
  long time;

  /**
   * Constructor.
   *
   * @param size The size of the DFT.
   */
  public TornadoDFT(int size) {
    inReal = new float[size];
    inImag = new float[size];
    outReal = new float[size];
    outImag = new float[size];
    for (int i = 0; i < size; i++) {
      inReal[i] = 1 / (float) (i + 2);
      inImag[i] = 1 / (float) (i + 2);
    }
  }

  /**
   * Compute the Discrete Fourier Transform.
   *
   * @param inreal  Input real values.
   * @param inimag  Input imaginary values.
   * @param outreal Output real values.
   * @param outimag Output imaginary values.
   */
  public static void computeDft(float[] inreal, float[] inimag, float[] outreal, float[] outimag) {
    int n = inreal.length;
    for (@Parallel int k = 0; k < n; k++) { // For each output element
      float sumReal = 0;
      float simImag = 0;
      for (int t = 0; t < n; t++) { // For each input element
        float angle = (2 * floatPI() * t * k) / n;
        sumReal += inreal[t] * cos(angle) + inimag[t] * sin(angle);
        simImag += -inreal[t] * sin(angle) + inimag[t] * cos(angle);
      }
      outreal[k] = sumReal;
      outimag[k] = simImag;
    }
  }

  /**
   * Execute the Discrete Fourier Transform on a TornadoDevice.
   *
   * @param device The TornadoDevice to use.
   */
  public void execute(TornadoDevice device) {
    TaskGraph graph =
        new TaskGraph("DFT").transferToDevice(EVERY_EXECUTION, inReal, inImag)
            .task("t0", TornadoDFT::computeDft, inReal, inImag, outReal, outImag)
            .transferToHost(EVERY_EXECUTION, outReal, outImag);

    ImmutableTaskGraph itg = graph.snapshot();
    TornadoExecutionPlan executionPlan = new TornadoExecutionPlan(itg);
    executionPlan.withWarmUp().withDevice(device);
    time = -System.nanoTime();
    executionPlan.execute();
    time += System.nanoTime();
  }

  /**
   * Execute the Discrete Fourier Transform on the default TornadoDevice.
   */
  public void execute() {
    TornadoDevice device = TornadoRuntimeProvider.getTornadoRuntime().getDefaultDevice();
    execute(device);
  }

  /**
   * Validate the Discrete Fourier Transform on the default TornadoDevice.
   *
   * @param deviceID The device ID to use.
   */
  public void validate(int deviceID) {
    TornadoDevice device = FFXTornado.getDevice(deviceID);
    validate(device);
  }

  /**
   * Validate the Discrete Fourier Transform on a TornadoDevice.
   *
   * @param device The TornadoDevice to use.
   */
  public void validate(TornadoDevice device) {
    execute(device);

    long javaTime = -System.nanoTime();
    computeDft(inReal, inImag, outReal, outImag);
    javaTime += System.nanoTime();

    System.out.println(" ");
    FFXTornado.logDevice(device);
    double speedUp = (double) javaTime / (double) time;
    System.out.printf(" %12s %8.6f (sec)\n %12s %8.6f (sec) Speed-Up %8.6f%n",
        " Java", 1.0e-9 * javaTime, " OpenCL", 1.0e-9 * time, speedUp);

  }
}
