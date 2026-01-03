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
package ffx.openmm;

import com.sun.jna.Memory;
import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_findPlatform;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getDefaultPluginsDirectory;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getNumPlatforms;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getOpenMMVersion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPlatform;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPlatformByName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPlatform_1;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPluginLoadFailures;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPropertyDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPropertyNames;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPropertyValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getSpeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_loadPluginLibrary;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_loadPluginsFromDirectory;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_registerPlatform;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_setPropertyDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_setPropertyValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_supportsDoublePrecision;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_supportsKernels;

/**
 * A Platform defines an implementation of all the kernels needed to perform some calculation.
 * More precisely, a Platform object acts as a registry for a set of KernelFactory
 * objects which together implement the kernels.  The Platform class, in turn, provides a
 * static registry of all available Platform objects.
 */
public class Platform {

  /**
   * OpenMM Platform pointer.
   */
  private PointerByReference pointer;

  /**
   * OpenMM Platform constructor.
   *
   * @param pointer The OpenMM Platform pointer.
   */
  public Platform(PointerByReference pointer) {
    this.pointer = pointer;
  }

  /**
   * OpenMM Platform constructor.
   *
   * @param platformName The name of the OpenMM Platform.
   */
  public Platform(String platformName) {
    pointer = OpenMM_Platform_getPlatformByName(platformName);
  }

  /**
   * Default constructor.
   */
  public Platform() {
    pointer = null;
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
   * Find a platform that supports a specified set of kernels.
   *
   * @param kernelNames The set of kernels that must be supported.
   * @return A platform that supports the specified kernels.
   */
  public static Platform findPlatform(StringArray kernelNames) {
    PointerByReference platformPointer = OpenMM_Platform_findPlatform(kernelNames.getPointer());
    return new Platform(platformPointer);
  }

  /**
   * Get the default directory from which to load plugins.
   *
   * @return The default directory from which to load plugins.
   */
  public static String getDefaultPluginsDirectory() {
    Pointer directory = OpenMM_Platform_getDefaultPluginsDirectory();
    return directory.getString(0);
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
   * Get the OpenMM Platform pointer.
   *
   * @return The OpenMM Platform pointer.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

  /**
   * Get a registered platform by index.
   *
   * @param index The index of the platform to get.
   * @return The platform with the specified index.
   */
  public static Platform getPlatform(int index) {
    PointerByReference platformPointer = OpenMM_Platform_getPlatform(index);
    return new Platform(platformPointer);
  }

  /**
   * Get a registered platform by name.
   *
   * @param name The name of the platform to get.
   * @return The platform with the specified name.
   */
  public static Platform getPlatform_1(String name) {
    PointerByReference platformPointer = OpenMM_Platform_getPlatform_1(name);
    return new Platform(platformPointer);
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
   * Get the default value of a platform property.
   *
   * @param property The name of the property to get.
   * @return The default value of the property.
   */
  public String getPropertyDefaultValue(String property) {
    Pointer propertyPointer = pointerForString(property);
    Pointer value = OpenMM_Platform_getPropertyDefaultValue(pointer, propertyPointer);
    if (value == null) {
      return null;
    }
    return value.getString(0);
  }

  /**
   * Get the names of all platform properties.
   *
   * @return A StringArray containing the names of all platform properties.
   */
  public StringArray getPropertyNames() {
    return new StringArray(OpenMM_Platform_getPropertyNames(pointer));
  }

  /**
   * Get the value of a context-specific platform property.
   *
   * @param context  The context for which to get the property value.
   * @param property The name of the property to get.
   * @return The value of the property.
   */
  public String getPropertyValue(Context context, String property) {
    Pointer propertyPointer = pointerForString(property);
    Pointer value = OpenMM_Platform_getPropertyValue(pointer, context.getPointer(), propertyPointer);
    if (value == null) {
      return null;
    }
    return value.getString(0);
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
   * Load plugins from a directory.
   *
   * @param directory The directory to load plugins from.
   * @return The OpenMMStringArray of plugins loaded.
   */
  public static StringArray loadPluginsFromDirectory(String directory) {
    return new StringArray(OpenMM_Platform_loadPluginsFromDirectory(directory));
  }

  /**
   * Load a dynamic library that contains a plugin.
   *
   * @param file The path to the dynamic library file.
   */
  public static void loadPluginLibrary(String file) {
    OpenMM_Platform_loadPluginLibrary(file);
  }

  /**
   * Register a new platform.
   *
   * @param platform The platform to register.
   */
  public static void registerPlatform(Platform platform) {
    OpenMM_Platform_registerPlatform(platform.getPointer());
  }

  /**
   * Set the default value of a platform property.
   *
   * @param property The name of the property to set.
   * @param value    The new default value for the property.
   */
  public void setPropertyDefaultValue(String property, String value) {
    Pointer propertyPointer = pointerForString(property);
    Pointer valuePointer = pointerForString(value);
    OpenMM_Platform_setPropertyDefaultValue(pointer, propertyPointer, valuePointer);
  }

  /**
   * Set the value of a context-specific platform property.
   *
   * @param context  The context for which to set the property value.
   * @param property The name of the property to set.
   * @param value    The new value for the property.
   */
  public void setPropertyValue(Context context, String property, String value) {
    Pointer propertyPointer = pointerForString(property);
    Pointer valuePointer = pointerForString(value);
    OpenMM_Platform_setPropertyValue(pointer, context.getPointer(), propertyPointer, valuePointer);
  }

  /**
   * Determine whether this Platform supports double precision arithmetic.
   *
   * @return true if the Platform supports double precision, false otherwise.
   */
  public boolean supportsDoublePrecision() {
    int result = OpenMM_Platform_supportsDoublePrecision(pointer);
    return result != 0;
  }

  /**
   * Determine whether this Platform supports a specified set of kernels.
   *
   * @param kernelNames The set of kernels to test for support.
   * @return true if the Platform supports all of the specified kernels, false otherwise.
   */
  public boolean supportsKernels(StringArray kernelNames) {
    int result = OpenMM_Platform_supportsKernels(pointer, kernelNames.getPointer());
    return result != 0;
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