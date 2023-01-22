// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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

import static java.lang.String.format;
import static uk.ac.manchester.tornado.api.collections.math.TornadoMath.cos;
import static uk.ac.manchester.tornado.api.collections.math.TornadoMath.floatPI;
import static uk.ac.manchester.tornado.api.collections.math.TornadoMath.sin;

import ffx.numerics.tornado.FFXTornado;
import java.util.logging.Logger;
import uk.ac.manchester.tornado.api.TaskSchedule;
import uk.ac.manchester.tornado.api.annotations.Parallel;
import uk.ac.manchester.tornado.api.common.TornadoDevice;
import uk.ac.manchester.tornado.api.runtime.TornadoRuntime;

/** Proof-of-concept use of the TornadoVM for parallelization of Java code. */
public class TornadoDFT {

  private static final Logger logger = Logger.getLogger(TornadoDFT.class.getName());
  float[] inReal;
  float[] inImag;
  float[] outReal;
  float[] outImag;
  long time;
  private int size;

  public TornadoDFT(int size) {
    this.size = size;
    inReal = new float[size];
    inImag = new float[size];
    outReal = new float[size];
    outImag = new float[size];
    for (int i = 0; i < size; i++) {
      inReal[i] = 1 / (float) (i + 2);
      inImag[i] = 1 / (float) (i + 2);
    }
  }

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

  public void execute(TornadoDevice device) {
    TaskSchedule graph =
        new TaskSchedule("DFT")
            .streamIn(inReal, inImag)
            .task("t0", TornadoDFT::computeDft, inReal, inImag, outReal, outImag)
            .streamOut(outReal, outImag);

    graph.setDevice(device);
    graph.warmup();

    time = -System.nanoTime();
    graph.execute();
    time += System.nanoTime();
  }

  public void execute() {
    TornadoDevice device = TornadoRuntime.getTornadoRuntime().getDefaultDevice();
    execute(device);
  }

  public void validate(int deviceID) {
    TornadoDevice device = FFXTornado.getDevice(deviceID);
    validate(device);
  }

  public void validate(TornadoDevice device) {
    execute(device);

    long javaTime = -System.nanoTime();
    computeDft(inReal, inImag, outReal, outImag);
    javaTime += System.nanoTime();

    System.out.println(" ");
    FFXTornado.logDevice(device);
    double speedUp = (double) javaTime / (double) time;
    System.out.println(format(" %12s %8.6f (sec)\n %12s %8.6f (sec) Speed-Up %8.6f",
        " Java", 1.0e-9 * javaTime, " OpenCL", 1.0e-9 * time, speedUp));

  }
}
