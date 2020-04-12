//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
//******************************************************************************
package ffx.utilities;

import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;

/**
 * <p>Abstract BaseFFXTest class.</p>
 *
 * @author Michael J. Schnieders
 */
public abstract class BaseFFXTest {
    public static final boolean ffxCI = System.getProperty("ffx.ci", "false").equalsIgnoreCase("true");
    public static final boolean ffxOpenMM = System.getProperty("ffx.openMM", "false").equalsIgnoreCase("true");
    /**
     * Constant <code>logger</code>
     */
    protected static final Logger logger = Logger.getLogger(BaseFFXTest.class.getName());
    private static final Level origLevel = Logger.getLogger("ffx").getLevel();
    private static final Level testLevel;
    private static final Level ffxLevel;
    private static Properties properties;

    static {
        Level level;
        try {
            level = Level.parse(System.getProperty("ffx.test.log", "INFO").toUpperCase());
        } catch (Exception ex) {
            logger.warning(String.format(" Exception %s in parsing value of ffx.test.log", ex));
            level = origLevel;
        }
        testLevel = level;

        try {
            level = Level.parse(System.getProperty("ffx.log", "INFO").toUpperCase());
        } catch (Exception ex) {
            logger.warning(String.format(" Exception %s in parsing value of ffx.log", ex));
            level = origLevel;
        }
        ffxLevel = level;
    }

    /**
     * <p>afterClass.</p>
     */
    @AfterClass
    public static void afterClass() {
        Logger.getLogger("ffx").setLevel(origLevel);
        logger.setLevel(origLevel);
    }

    /**
     * <p>afterTest.</p>
     */
    @After
    public void afterTest() {
        // All properties are set to the values they were at the beginning of the test.
        System.setProperties(properties);
    }

    /**
     * <p>beforeClass.</p>
     */
    @BeforeClass
    public static void beforeClass() {
        // Set appropriate logging levels for interior/exterior Loggers.
        Logger.getLogger("ffx").setLevel(ffxLevel);
        logger.setLevel(testLevel);
    }

    /**
     * <p>beforeTest.</p>
     */
    @Before
    public void beforeTest() {
        // New properties object that will hold the property key-value pairs that were present
        // at the beginning of the test.
        properties = new Properties();

        // currentProperties holds the properties at the beginning of the test.
        Properties currentProperties = System.getProperties();

        // All key-value pairs from currentProperties are stored in the properties object.
        for (String key : currentProperties.stringPropertyNames()) {
            properties.setProperty(key, currentProperties.getProperty(key));
        }
    }
}
