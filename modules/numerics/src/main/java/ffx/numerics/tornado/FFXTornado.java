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
package ffx.numerics.tornado;

import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Logger;
import uk.ac.manchester.tornado.api.TornadoDeviceContext;
import uk.ac.manchester.tornado.api.TornadoDriver;
import uk.ac.manchester.tornado.api.TornadoTargetDevice;
import uk.ac.manchester.tornado.api.common.Event;
import uk.ac.manchester.tornado.api.common.TornadoDevice;
import uk.ac.manchester.tornado.api.enums.TornadoDeviceType;
import uk.ac.manchester.tornado.api.enums.TornadoVMBackendType;
import uk.ac.manchester.tornado.api.memory.TornadoDeviceObjectState;
import uk.ac.manchester.tornado.api.memory.TornadoMemoryProvider;
import uk.ac.manchester.tornado.api.runtime.TornadoRuntime;

/** Utility Routines to use the TornadoVM */
public class FFXTornado {

  private static final Logger logger = Logger.getLogger(FFXTornado.class.getName());

  /**
   * Get the default Tornado Device.
   *
   * @return The default TornadoDevice instance.
   */
  public static TornadoDevice getDevice() {
    return TornadoRuntime.getTornadoRuntime().getDefaultDevice();
  }

  /**
   * Get the Tornado Device using specified driver and device index.
   *
   * @param driverIndex Driver index.
   * @param deviceIndex Device index.
   * @return The TornadoDevice instance.
   */
  public static TornadoDevice getDevice(int driverIndex, int deviceIndex) {
    TornadoDriver tornadoDriver = TornadoRuntime.getTornadoRuntime().getDriver(driverIndex);
    return tornadoDriver.getDevice(deviceIndex);
  }

  /**
   * Get the specified Tornado Device.
   *
   * @param deviceID The device ID.
   * @return The TornadoDevice instance.
   */
  public static TornadoDevice getDevice(int deviceID) {
    int n = 0;
    int numDrivers = TornadoRuntime.getTornadoRuntime().getNumDrivers();
    for (int driverIndex = 0; driverIndex < numDrivers; driverIndex++) {
      TornadoDriver driver = TornadoRuntime.getTornadoRuntime().getDriver(driverIndex);
      for (int deviceIndex = 0; deviceIndex < driver.getDeviceCount(); deviceIndex++) {
        if (n == deviceID) {
          TornadoRuntime.setProperty("devices", driverIndex + ":" + deviceIndex);
          return getDevice(driverIndex, deviceIndex);
        }
        n++;
      }
    }
    return null;
  }

  /**
   * Get all TornadoDevice instances.
   *
   * @return A List of TornadoDevice instances.
   */
  public static int getNumberOfDevices() {
    int n = 0;
    int numDrivers = TornadoRuntime.getTornadoRuntime().getNumDrivers();
    for (int driverIndex = 0; driverIndex < numDrivers; driverIndex++) {
      TornadoDriver driver = TornadoRuntime.getTornadoRuntime().getDriver(driverIndex);
      int count = driver.getDeviceCount();
      n += count;
    }
    return n;
  }

  /**
   * List details about the passed TornadoDevice instance.
   *
   * @param device The TornadoDevice instance.
   */
  public static void logDevice(TornadoDevice device) {
    TornadoTargetDevice tornadoTargetDevice = device.getPhysicalDevice();
    long[] workItemSize = tornadoTargetDevice.getDeviceMaxWorkItemSizes();
    //        logger.info(format("\n Device Name:         %s",
    // tornadoTargetDevice.getDeviceName()));
    //        logger.info(format(" Compute Units:       %s",
    // tornadoTargetDevice.getDeviceMaxComputeUnits()));
    //        logger.info(format(" Max Work Item Sizes: [%d, %d, %d]", workItemSize[0],
    // workItemSize[1], workItemSize[2]));
    //        logger.info(format(" Clock Frequency:     %6d Ghz",
    // tornadoTargetDevice.getDeviceMaxClockFrequency()));
    //        logger.info(format(" Global Memory:       %6d MB",
    // tornadoTargetDevice.getDeviceGlobalMemorySize() / 1024 / 1024));
    //        logger.info(format(" Local Memory:        %6d KB",
    // tornadoTargetDevice.getDeviceLocalMemorySize() / 1024));
    System.out.printf("\n Device Name:         %s%n", tornadoTargetDevice.getDeviceName());
    System.out.printf(" Backend:             %s%n", device.getTornadoVMBackend().name());
    System.out.printf(" Compute Units:       %s%n", tornadoTargetDevice.getDeviceMaxComputeUnits());
    System.out.printf(" Max Work Item Sizes: [%d, %d, %d]%n", workItemSize[0], workItemSize[1],
        workItemSize[2]);
    System.out.printf(" Clock Frequency:     %6d Ghz%n",
        tornadoTargetDevice.getDeviceMaxClockFrequency());
    System.out.printf(" Global Memory:       %6d MB%n",
        tornadoTargetDevice.getDeviceGlobalMemorySize() / 1024 / 1024);
    System.out.printf(" Local Memory:        %6d KB%n",
        tornadoTargetDevice.getDeviceLocalMemorySize() / 1024);
  }
}
