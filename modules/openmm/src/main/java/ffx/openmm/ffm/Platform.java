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
package ffx.openmm.ffm;

import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.MemorySegment;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * FFM-backed OpenMM platform, which selects and configures computational kernels.
 *
 * <p>OpenMM's platform registry owns the platform objects. This wrapper borrows a registry handle,
 * so {@link #destroy()} invalidates only this Java wrapper and never destroys a native platform.</p>
 */
public class Platform extends OpenMMHandle {

  /**
   * Create an uninitialized platform wrapper for deferred setup.
   *
   * <p>Native operations require a platform handle; this constructor exists for compatibility with
   * the JNA façade.</p>
   */
  public Platform() {
    super(MemorySegment.NULL, true);
  }

  /**
   * @param name name of a registered platform; the native registry lookup must succeed.
   */
  public Platform(String name) {
    super(byName(name));
  }

  /**
   * Wrap a registered native platform handle without taking ownership of the platform.
   *
   * @param pointer borrowed native platform handle.
   */
  public Platform(MemorySegment pointer) {
    super(pointer);
  }

  /**
   * Invalidate this borrowed platform wrapper without destroying the registered native platform.
   */
  @Override
  public void destroy() {
    invalidate();
  }

  /**
   * Find the fastest registered platform that supports all requested kernels.
   *
   * <p>OpenMM reports a failure when no registered platform supports every requested kernel.</p>
   *
   * @param kernelNames non-null array of required kernel names; ownership remains with the caller.
   * @return borrowed wrapper for the fastest compatible registered platform.
   */
  public static Platform findPlatform(StringArray kernelNames) {
    Objects.requireNonNull(kernelNames, "Kernel names cannot be null.");
    initialize();
    return new Platform(OpenMMNative.OpenMM_Platform_findPlatform(kernelNames.getPointer()));
  }

  /**
   * Get OpenMM's default plugin directory, honoring OpenMM's own {@code OPENMM_PLUGIN_DIR}
   * configuration and platform-specific default.
   *
   * @return copied plugin-directory path.
   */
  public static String getDefaultPluginsDirectory() {
    initialize();
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_Platform_getDefaultPluginsDirectory.makeInvoker().apply());
  }

  /**
   * @return copied registered platform name.
   */
  public String getName() {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_Platform_getName(getPointer()));
  }

  /**
   * @return number of registered platforms.
   */
  public static int getNumPlatforms() {
    initialize();
    return OpenMMNative.OpenMM_Platform_getNumPlatforms.makeInvoker().apply();
  }

  /**
   * @return copied OpenMM version string.
   */
  public static String getOpenMMVersion() {
    initialize();
    return OpenMMStrings.copy(OpenMMNative.OpenMM_Platform_getOpenMMVersion.makeInvoker().apply());
  }

  /**
   * Get a registered platform by its zero-based registry index.
   *
   * @param index zero-based platform index.
   * @return borrowed wrapper for the registered platform.
   */
  public static Platform getPlatform(int index) {
    initialize();
    return new Platform(OpenMMNative.OpenMM_Platform_getPlatform(index));
  }

  /**
   * Get a registered platform by name.
   *
   * <p>OpenMM reports a failure if no platform with that name is registered.</p>
   *
   * @param name non-null registered platform name.
   * @return borrowed wrapper for the registered platform.
   */
  public static Platform getPlatform(String name) {
    return new Platform(byName(name));
  }

  /**
   * Legacy name-based lookup alias retained for source compatibility with the JNA façade.
   *
   * @param name registered platform name.
   * @return platform wrapper.
   * @deprecated use {@link #getPlatform(String)}.
   */
  @Deprecated
  public static Platform getPlatform_1(String name) {
    return getPlatform(name);
  }

  /**
   * Copy and release the native string array of failures from OpenMM's most recent plugin-load
   * operation.
   *
   * @return copied failure messages.
   */
  public static String[] getPluginLoadFailures() {
    initialize();
    return copyAndDestroy(
        OpenMMNative.OpenMM_Platform_getPluginLoadFailures.makeInvoker().apply());
  }

  /**
   * Get the default value of a platform-specific property.
   *
   * @param property non-null property name.
   * @return copied property value.
   */
  public String getPropertyDefaultValue(String property) {
    return OpenMMStrings.withUtf8StringResult(property,
        value -> OpenMMStrings.copy(
            OpenMMNative.OpenMM_Platform_getPropertyDefaultValue(getPointer(), value)));
  }

  /**
   * Get the names of platform-specific properties supported by this platform.
   *
   * @return copied property names.
   */
  public String[] getPropertyNames() {
    return copyAndDestroy(OpenMMNative.OpenMM_Platform_getPropertyNames(getPointer()));
  }

  /**
   * Get the value of a platform-specific property for a context.
   *
   * @param context live context whose platform property is queried.
   * @param property non-null property name.
   * @return copied property value.
   */
  public String getPropertyValue(Context context, String property) {
    Objects.requireNonNull(context, "Context cannot be null.");
    return OpenMMStrings.withUtf8StringResult(property,
        value -> OpenMMStrings.copy(
            OpenMMNative.OpenMM_Platform_getPropertyValue(
                getPointer(), context.getPointer(), value)));
  }

  /**
   * @return relative platform speed estimate; OpenMM's reference implementation uses 1.0 and
   *         other platforms report an estimate relative to it.
   */
  public double getSpeed() {
    return OpenMMNative.OpenMM_Platform_getSpeed(getPointer());
  }

  /**
   * Load plugins found in a platform-delimited list of directories.
   *
   * <p>On Unix-like systems the native API uses {@code :} between paths; on Windows it uses
   * {@code ;}. This method calls OpenMM directly and does not apply the RPMD-plugin filtering used
   * by {@link OpenMMRuntime#initialize()}. Initialize the OpenMM runtime before calling this
   * method; it does not load the base libraries itself.</p>
   *
   * @param directory non-null platform-delimited plugin-directory list.
   * @return copied names of libraries successfully loaded by OpenMM.
   */
  public static String[] loadPluginsFromDirectory(String directory) {
    return OpenMMStrings.withUtf8StringResult(directory,
        value -> copyAndDestroy(OpenMMNative.OpenMM_Platform_loadPluginsFromDirectory(value)));
  }

  /**
   * Load a plugin shared library through OpenMM.
   *
   * <p>Initialize the OpenMM runtime before calling this method; it does not load the base
   * libraries itself.</p>
   *
   * @param file non-null plugin shared-library path accepted by the platform's native loader.
   */
  public static void loadPluginLibrary(String file) {
    OpenMMStrings.withUtf8String(file, OpenMMNative::OpenMM_Platform_loadPluginLibrary);
  }

  /**
   * Register a platform with OpenMM's process-wide registry.
   *
   * @param platform live platform wrapper to register.
   */
  public static void registerPlatform(Platform platform) {
    Objects.requireNonNull(platform, "Platform cannot be null.");
    initialize();
    OpenMMNative.OpenMM_Platform_registerPlatform(platform.getPointer());
  }

  /**
   * Set the default value of a platform-specific property for newly created contexts.
   *
   * @param property non-null property name.
   * @param value    non-null value to use as the default.
   */
  public void setPropertyDefaultValue(String property, String value) {
    OpenMMStrings.withUtf8String(property, propertyPointer ->
        OpenMMStrings.withUtf8String(value, valuePointer ->
            OpenMMNative.OpenMM_Platform_setPropertyDefaultValue(
                getPointer(), propertyPointer, valuePointer)));
  }

  /**
   * Set the value of a platform-specific property for a context.
   *
   * @param context live context to configure.
   * @param property non-null property name.
   * @param value    non-null property value.
   */
  public void setPropertyValue(Context context, String property, String value) {
    Objects.requireNonNull(context, "Context cannot be null.");
    OpenMMStrings.withUtf8String(property, propertyPointer ->
        OpenMMStrings.withUtf8String(value, valuePointer ->
            OpenMMNative.OpenMM_Platform_setPropertyValue(
                getPointer(), context.getPointer(), propertyPointer, valuePointer)));
  }

  /**
   * Report OpenMM's legacy double-precision support flag.
   *
   * <p>This query is deprecated by the OpenMM API because a single boolean does not describe
   * platforms with multiple precision modes.</p>
   *
   * @return whether this platform reports double-precision support.
   */
  public boolean supportsDoublePrecision() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_Platform_supportsDoublePrecision(getPointer()));
  }

  /**
   * Determine whether this platform implements every requested kernel.
   *
   * @param kernelNames non-null array of required kernel names; ownership remains with the caller.
   * @return {@code true} if all requested kernels are supported.
   */
  public boolean supportsKernels(StringArray kernelNames) {
    Objects.requireNonNull(kernelNames, "Kernel names cannot be null.");
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_Platform_supportsKernels(getPointer(), kernelNames.getPointer()));
  }

  /**
   * Resolve a registered platform by name.
   */
  private static MemorySegment byName(String name) {
    return OpenMMStrings.withUtf8StringResult(name, value -> {
      initialize();
      return OpenMMNative.OpenMM_Platform_getPlatformByName(value);
    });
  }

  /**
   * Copy and release an OpenMM-owned string array returned to the caller.
   */
  private static String[] copyAndDestroy(MemorySegment array) {
    try {
      int size = OpenMMNative.OpenMM_StringArray_getSize(array);
      List<String> values = new ArrayList<>(size);
      for (int index = 0; index < size; index++) {
        values.add(OpenMMStrings.copy(OpenMMNative.OpenMM_StringArray_get(array, index)));
      }
      return values.toArray(String[]::new);
    } finally {
      OpenMMNative.OpenMM_StringArray_destroy(array);
    }
  }

  /**
   * Initialize the FFM native runtime before resolving a platform.
   */
  private static void initialize() {
    OpenMMRuntime.initialize();
  }
}
