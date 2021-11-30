// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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

import static java.lang.String.format;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Properties;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;

/**
 * The FFXTest configures the context for FFX tests. This includes: 1) Sets testing related
 * environment variables:
 * <br>
 * Set "-Dffx.ci=true" for a CI environment (default: false).
 * <br>
 * Set "-Dffx.openMM=true" for CUDA dependent OpenMM tests (default: false).
 * <br>
 * 2) Configures the logging level, using the ffx.test.log System property (default: INFO).
 * <br>
 * 3) Stores System properties prior to each test, and restores them after each test (i.e. properties
 * set by a test do not effect the next test).
 * <br>
 * 4) Upon request, creates a temporary that is deleted after the current test is completed.
 * <br>
 * 5) Requests garbage collection after each test.
 *
 * @author Michael J. Schnieders
 */
public abstract class FFXTest {

  /** Constant <code>ffxCI=System.getProperty("ffx.ci", "false").equalsIgnoreCase("true")</code> */
  public static final boolean ffxCI =
      System.getProperty("ffx.ci", "false").equalsIgnoreCase("true");

  /**
   * Constant <code>ffxOpenMM=System.getProperty("ffx.openMM", "false").equalsIgnoreCase("true")
   * </code>
   */
  public static final boolean ffxOpenMM =
      System.getProperty("ffx.openMM", "false").equalsIgnoreCase("true");

  /** Constant <code>logger</code> */
  protected static final Logger logger = Logger.getLogger(FFXTest.class.getName());

  private static final Level origLevel;
  private static final Level testLevel;
  private static Properties properties;
  private Path path = null;

  static {
    Level level;
    try {
      level = Level.parse(System.getProperty("ffx.log", "INFO").toUpperCase());
    } catch (Exception ex) {
      logger.warning(format(" Exception %s in parsing value of ffx.log", ex));
      level = Level.INFO;
    }
    origLevel = level;

    try {
      level = Level.parse(System.getProperty("ffx.test.log", "WARNING").toUpperCase());
    } catch (Exception ex) {
      logger.warning(format(" Exception %s in parsing value of ffx.test.log", ex));
      level = origLevel;
    }
    testLevel = level;
    System.setProperty("ffx.log", testLevel.toString());
  }

  /** afterClass. */
  @AfterClass
  public static void afterClass() {
    Logger.getGlobal().setLevel(origLevel);
    Logger.getLogger("ffx").setLevel(origLevel);
    logger.setLevel(origLevel);
  }

  /** beforeClass. */
  @BeforeClass
  public static void beforeClass() {
    // Set appropriate logging levels for interior/exterior Loggers.
    Logger.getGlobal().setLevel(testLevel);
    Logger.getLogger("ffx").setLevel(testLevel);
    logger.setLevel(testLevel);
  }

  /**
   * Create temporary testing directory that will be deleted after the current test.
   *
   * @return Path to the testing directory.
   */
  public Path registerTemporaryDirectory() {
    deleteTemporaryDirectory();
    try {
      path = Files.createTempDirectory("FFXTestDirectory");
    } catch (java.io.IOException e) {
      fail(" Could not create a temporary directory.");
    }
    return path;
  }

  /**
   * Delete the temporary test directory.
   */
  private void deleteTemporaryDirectory() {
    // Delete the test directory if it exists.
    if (path != null) {
      try {
        DirectoryUtils.deleteDirectoryTree(path);
      } catch (IOException e) {
        System.out.println(e.toString());
        fail(" Exception deleting files created by Frac2Cart.");
      }
      path = null;
    }
  }

  /** afterTest. */
  @After
  public void afterTest() {
    // All properties are set to the values they were at the beginning of the test.
    System.setProperties(properties);

    // Delete the test directory if it exists.
    deleteTemporaryDirectory();

    // Collect garbage.
    System.gc();
  }

  /** beforeTest. */
  @Before
  public void beforeTest() {
    // New properties object that will hold the property key-value pairs that were present
    // at the beginning of the test.
    properties = new Properties();

    // currentProperties holds the properties at the beginning of the test.
    Properties currentProperties = System.getProperties();

    // All key-value pairs from currentProperties are stored in the properties object.
    currentProperties.stringPropertyNames()
        .forEach(key -> properties.setProperty(key, currentProperties.getProperty(key)));
  }
}
