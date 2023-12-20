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
package ffx.openmm;

import com.sun.jna.Memory;
import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;
import edu.uiowa.jopenmm.OpenMMLibrary;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getNumPlatforms;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getOpenMMVersion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPlatformByName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPluginLoadFailures;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getSpeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_loadPluginsFromDirectory;

/**
 * OpenMM Platform.
 */
public class Platform {

  /**
   * OpenMM Platform pointer.
   */
  private PointerByReference pointer;

  /**
   * OpenMM Platform constructor.
   *
   * @param platformName The name of the OpenMM Platform.
   */
  public Platform(String platformName) {
    pointer = OpenMM_Platform_getPlatformByName(platformName);
  }

  /**
   * Get the OpenMM Platform pointer.
   *
   * @return The OpenMM Platform pointer.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

  /**
   * Get the name of the OpenMM Platform.
   *
   * @return The name of the OpenMM Platform.
   */
  public String getName() {
    Pointer name = OpenMM_Platform_getName(pointer);
    return name.getString(0);
  }

  /**
   * Get an estimate of how fast this Platform class is.
   *
   * @return The speed of the OpenMM Platform.
   */
  public double getSpeed() {
    return OpenMM_Platform_getSpeed(pointer);
  }

  /**
   * Set an OpenMM Platform property.
   */
  public void setPropertyDefaultValue(String propertyName, String defaultValue) {
    Pointer name = pointerForString(propertyName);
    Pointer value = pointerForString(defaultValue);
    OpenMMLibrary.OpenMM_Platform_setPropertyDefaultValue(pointer, name, value);
  }

  /**
   * Get the number of OpenMM Platforms.
   *
   * @return The number of OpenMM Platforms.
   */
  public static int getNumPlatforms() {
    return OpenMM_Platform_getNumPlatforms();
  }

  /**
   * Get the OpenMM version.
   *
   * @return The version of OpenMM.
   */
  public static String getOpenMMVersion() {
    Pointer version = OpenMM_Platform_getOpenMMVersion();
    return version.getString(0);
  }

  /**
   * Load plugins from a directory.
   *
   * @param directory The directory to load plugins from.
   * @return The OpenMMStringArray of plugins loaded.
   */
  public static StringArray loadPluginsFromDirectory(String directory) {
    return new StringArray(OpenMM_Platform_loadPluginsFromDirectory(directory));
  }

  /**
   * Get the plugin load failures.
   *
   * @return The OpenMMStringArray of plugin load failures.
   */
  public static StringArray getPluginLoadFailures() {
    return new StringArray(OpenMM_Platform_getPluginLoadFailures());
  }

  /**
   * Destroy the OpenMM Platform instance.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_Platform_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Create a JNA Pointer to a String.
   *
   * @param string WARNING: assumes ascii-only string
   * @return pointer.
   */
  private static Pointer pointerForString(String string) {
    Pointer pointer = new Memory(string.length() + 1);
    pointer.setString(0, string);
    return pointer;
  }

}
