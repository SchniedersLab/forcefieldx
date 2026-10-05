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

import java.io.IOException;
import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

/**
 * Loads the OpenMM native libraries and supported plugins for the generated FFM bindings.
 *
 * <p>Initialization requires {@value #LIBRARY_DIRECTORY_ENVIRONMENT} to identify a directory
 * containing the platform-mapped OpenMM, OpenMMAmoeba, and OpenMMDrude libraries. Plugins are read
 * from {@value #PLUGIN_DIRECTORY_ENVIRONMENT}, or from the {@code plugins} subdirectory of the
 * library directory when that optional variable is unset or blank. RPMD plugin files are excluded
 * because this FFM runtime does not support the RPMD extension.</p>
 *
 * <p>Initialization is process-wide and has no unload/reset operation. A plugin load failure
 * aborts initialization with an {@link IllegalStateException}.</p>
 */
public final class OpenMMRuntime {

  /**
   * Environment variable naming the directory containing the platform-mapped OpenMM, OpenMMAmoeba,
   * and OpenMMDrude shared libraries.
   */
  public static final String LIBRARY_DIRECTORY_ENVIRONMENT = "FFX_OPENMM_LIB_DIR";

  /**
   * Optional environment variable naming the OpenMM plugin directory. If unset or blank, the
   * {@code plugins} subdirectory of {@link #LIBRARY_DIRECTORY_ENVIRONMENT} is used.
   */
  public static final String PLUGIN_DIRECTORY_ENVIRONMENT = "FFX_OPENMM_PLUGIN_DIR";

  private static Path libraryDirectory;
  private static Path pluginDirectory;
  private static boolean initialized;

  /**
   * Prevent construction of this utility class.
   */
  private OpenMMRuntime() {
  }

  /**
   * Load the required OpenMM, AMOEBA, and Drude libraries and then load supported plugins.
   *
   * <p>This operation is idempotent after successful completion. It resolves configured
   * directories to real paths, excludes RPMD plugin files, and fails rather than completing when a
   * required library or plugin cannot be loaded.</p>
   *
   * @throws IllegalStateException if required directories or library files are missing, cannot be
   *                                resolved, or plugin loading reports failures.
   * @throws UnsatisfiedLinkError   if a required native library cannot be loaded.
   */
  public static synchronized void initialize() {
    if (initialized) {
      return;
    }

    Path configuredLibraryDirectory = requiredDirectory(
        LIBRARY_DIRECTORY_ENVIRONMENT, java.lang.System.getenv(LIBRARY_DIRECTORY_ENVIRONMENT));
    loadLibrary(configuredLibraryDirectory, "OpenMM");
    loadLibrary(configuredLibraryDirectory, "OpenMMAmoeba");
    loadLibrary(configuredLibraryDirectory, "OpenMMDrude");

    Path configuredPluginDirectory = pluginDirectory(configuredLibraryDirectory);
    Path filteredPluginDirectory = withoutRpmdPlugins(configuredPluginDirectory);
    loadPlugins(filteredPluginDirectory);

    libraryDirectory = configuredLibraryDirectory;
    pluginDirectory = filteredPluginDirectory;
    initialized = true;
  }

  /**
   * Return whether initialization has completed successfully.
   *
   * @return {@code true} after required libraries and the filtered plugin set have loaded.
   */
  public static synchronized boolean isInitialized() {
    return initialized;
  }

  /**
   * Get the resolved directory from which the OpenMM libraries were loaded.
   *
   * @return canonical native library directory.
   * @throws IllegalStateException if {@link #initialize()} has not completed successfully.
   */
  public static synchronized Path getLibraryDirectory() {
    requireInitialized();
    return libraryDirectory;
  }

  /**
   * Get the resolved plugin directory used by this runtime.
   *
   * <p>If RPMD plugin files were present, this is a temporary filtered directory containing links
   * to supported plugin files rather than the originally configured directory. The temporary
   * directory and links are scheduled for deletion when the process exits.</p>
   *
   * @return plugin directory passed to OpenMM's plugin loader.
   * @throws IllegalStateException if {@link #initialize()} has not completed successfully.
   */
  public static synchronized Path getPluginDirectory() {
    requireInitialized();
    return pluginDirectory;
  }

  /**
   * Load a required platform-mapped OpenMM library from a configured directory.
   */
  private static void loadLibrary(Path directory, String libraryName) {
    Path library = directory.resolve(java.lang.System.mapLibraryName(libraryName));
    if (!Files.isRegularFile(library)) {
      throw new IllegalStateException("Required OpenMM library does not exist: " + library);
    }
    java.lang.System.load(library.toString());
  }

  /**
   * Load plugins, failing explicitly when OpenMM reports a plugin-load failure.
   */
  private static void loadPlugins(Path directory) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment directoryString = arena.allocateFrom(directory.toString());
      readAndDestroy(OpenMMNative.OpenMM_Platform_loadPluginsFromDirectory(directoryString));
    }

    List<String> failures = readAndDestroy(
        OpenMMNative.OpenMM_Platform_getPluginLoadFailures.makeInvoker().apply());
    if (!failures.isEmpty()) {
      throw new IllegalStateException(
          "OpenMM plugins failed to load from " + directory + ": " + String.join("; ", failures));
    }
  }

  /**
   * Resolve the configured or default plugin directory.
   */
  private static Path pluginDirectory(Path configuredLibraryDirectory) {
    String configuredPluginDirectory = java.lang.System.getenv(PLUGIN_DIRECTORY_ENVIRONMENT);
    if (configuredPluginDirectory == null || configuredPluginDirectory.isBlank()) {
      return requiredDirectory(
          PLUGIN_DIRECTORY_ENVIRONMENT + " (default)", configuredLibraryDirectory.resolve("plugins").toString());
    }
    return requiredDirectory(PLUGIN_DIRECTORY_ENVIRONMENT, configuredPluginDirectory);
  }

  /**
   * Create a linked plugin directory without unsupported RPMD plugins when necessary.
   */
  private static Path withoutRpmdPlugins(Path configuredPluginDirectory) {
    List<Path> plugins;
    try (Stream<Path> paths = Files.list(configuredPluginDirectory)) {
      plugins = paths.filter(Files::isRegularFile).toList();
    } catch (IOException exception) {
      throw new IllegalStateException(
          "Cannot list OpenMM plugins in " + configuredPluginDirectory, exception);
    }

    if (plugins.stream().noneMatch(OpenMMRuntime::isRpmdPlugin)) {
      return configuredPluginDirectory;
    }

    try {
      Path filteredDirectory = Files.createTempDirectory("ffx-openmm-plugins-");
      filteredDirectory.toFile().deleteOnExit();
      for (Path plugin : plugins) {
        if (!isRpmdPlugin(plugin)) {
          Path link = filteredDirectory.resolve(plugin.getFileName());
          Files.createSymbolicLink(link, plugin);
          link.toFile().deleteOnExit();
        }
      }
      return filteredDirectory;
    } catch (IOException exception) {
      throw new IllegalStateException(
          "Cannot create an OpenMM plugin directory without RPMD plugins.", exception);
    }
  }

  /**
   * @return whether a plugin filename belongs to the unsupported RPMD extension.
   */
  private static boolean isRpmdPlugin(Path plugin) {
    return plugin.getFileName().toString().contains("OpenMMRPMD");
  }

  /**
   * Resolve and validate a required configuration directory.
   */
  private static Path requiredDirectory(String variableName, String configuredPath) {
    if (configuredPath == null || configuredPath.isBlank()) {
      throw new IllegalStateException("Set " + variableName + " to an OpenMM directory.");
    }

    try {
      Path directory = Path.of(configuredPath).toRealPath();
      if (!Files.isDirectory(directory)) {
        throw new IllegalStateException(variableName + " is not a directory: " + directory);
      }
      return directory;
    } catch (IOException exception) {
      throw new IllegalStateException(
          "Cannot resolve " + variableName + ": " + configuredPath, exception);
    }
  }

  /**
   * Copy and release an OpenMM-owned string array.
   */
  private static List<String> readAndDestroy(MemorySegment array) {
    try {
      int size = OpenMMNative.OpenMM_StringArray_getSize(array);
      List<String> values = new ArrayList<>(size);
      for (int index = 0; index < size; index++) {
        MemorySegment string = OpenMMNative.OpenMM_StringArray_get(array, index);
        values.add(string.reinterpret(Long.MAX_VALUE).getString(0));
      }
      return values;
    } finally {
      OpenMMNative.OpenMM_StringArray_destroy(array);
    }
  }

  /**
   * Ensure native library and plugin initialization has completed.
   */
  private static void requireInitialized() {
    if (!initialized) {
      throw new IllegalStateException("OpenMMRuntime.initialize() must be called first.");
    }
  }
}
