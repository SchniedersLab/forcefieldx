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
package ffx.utilities;

import static java.lang.Runtime.getRuntime;
import static java.lang.String.format;
import static java.lang.management.ManagementFactory.getOperatingSystemMXBean;

import com.sun.management.UnixOperatingSystemMXBean;
import java.lang.management.OperatingSystemMXBean;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Log resources.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Resources {

  private static final Logger logger = Logger.getLogger(Resources.class.getName());

  /**
   * Log resources.
   */
  public static void logResources() {
    if (logger.isLoggable(Level.INFO)) {
      StringBuilder sb = new StringBuilder("\n System Resources\n");

      Runtime runtime = getRuntime();
      runtime.runFinalization();
      runtime.gc();

      long MB = 1024 * 1024;
      OperatingSystemMXBean operatingSystemMXBean = getOperatingSystemMXBean();

      if (operatingSystemMXBean instanceof UnixOperatingSystemMXBean unixOS) {
        // CPU Time and Average System Load.
        long time = unixOS.getProcessCpuTime();
        sb.append(format("  Total CPU time:         %6.2f (sec)\n", time * 1.0e-9));
        double systemLoadAve = unixOS.getSystemLoadAverage();
        if (systemLoadAve >= 0) {
          sb.append(format("  System load average:    %6.2f\n", systemLoadAve));
        }
        // File handle use.
        long open = unixOS.getOpenFileDescriptorCount();
        long allowed = unixOS.getMaxFileDescriptorCount();
        sb.append(format("  Open file handles:      %6d of %6d allowed\n", open, allowed));

        // System Memory.
        long freePhysical = unixOS.getFreeMemorySize() / MB;
        long totalPhysical = unixOS.getTotalMemorySize() / MB;
        long freeSwap = unixOS.getFreeSwapSpaceSize() / MB;
        long totalSwap = unixOS.getTotalSwapSpaceSize() / MB;
        sb.append(format("  System memory:          %6d MB free out of %6d MB\n", freePhysical,
            totalPhysical));
        sb.append(
            format("  System swap space:      %6d MB free out of %6d MB\n", freeSwap, totalSwap));
      }

      // JVM Memory.
      sb.append(
          format("  JVM memory:             %6d MB free out of %6d MB", runtime.freeMemory() / MB,
              runtime.totalMemory() / MB));

      logger.info(sb.toString());
    }
  }
}
