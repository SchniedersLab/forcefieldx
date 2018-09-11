/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.utilities;

import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import com.sun.management.UnixOperatingSystemMXBean;

/**
 * Log resources.
 *
 * @author Michael J. Schnieders
 */
public class Resources {

    private static final Logger logger = Logger.getLogger(Resources.class.getName());

    /**
     * <p>logResources.</p>
     */
    public static void logResources() {
        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();

            Runtime runtime = Runtime.getRuntime();
            runtime.runFinalization();
            runtime.gc();

            long MB = 1024 * 1024;
            OperatingSystemMXBean os = ManagementFactory.getOperatingSystemMXBean();

            if (os instanceof UnixOperatingSystemMXBean) {
                UnixOperatingSystemMXBean unixOS = (UnixOperatingSystemMXBean) os;

                // CPU Time and Average System Load.
                Long time = unixOS.getProcessCpuTime();
                sb.append(format("\n Total CPU time:         %6.2f (sec)", time * 1.0e-9));
                double systemLoadAve = unixOS.getSystemLoadAverage();
                if (systemLoadAve >= 0) {
                    sb.append(format("\n System load average:    %6.2f", systemLoadAve));
                }
                // File handle use.
                long open = unixOS.getOpenFileDescriptorCount();
                long allowed = unixOS.getMaxFileDescriptorCount();
                sb.append(format("\n Open file handles:      %6d of %6d allowed",
                        open, allowed));

                // System Memory.
                long freePhysical = unixOS.getFreePhysicalMemorySize() / MB;
                long totalPhysical = unixOS.getTotalPhysicalMemorySize() / MB;
                long freeSwap = unixOS.getFreeSwapSpaceSize() / MB;
                long totalSwap = unixOS.getTotalSwapSpaceSize() / MB;
                sb.append(format("\n System memory:          %6d MB free out of %6d MB",
                        freePhysical, totalPhysical));
                sb.append(format("\n System swap space:      %6d MB free out of %6d MB",
                        freeSwap, totalSwap));

                // Log CPU usage.

                // double sysCPULoad = unixOS.getSystemCpuLoad() * 100.0;
                // double procCPULoad = unixOS.getProcessCpuLoad() * 100.0;
                /**
                 * if (procCPULoad >= 0) { sb.append(format("\n JVM CPU load:
                 * %6.3f%%", procCPULoad)); } if (sysCPULoad >= 0) {
                 * sb.append(format("\n System CPU load: %6.3f%%", sysCPULoad));
                 * }
                 */
            }

            // JVM Memory.
            sb.append(format("\n JVM memory:             %6d MB free out of %6d MB\n",
                    runtime.freeMemory() / MB, runtime.totalMemory() / MB));

            logger.info(sb.toString());
        }
    }
}
