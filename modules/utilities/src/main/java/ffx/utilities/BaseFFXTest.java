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
    /**
     * Constant <code>logger</code>
     */
    protected static final Logger logger = Logger.getLogger(BaseFFXTest.class.getName());
    private static final Level origLevel = Logger.getLogger("ffx").getLevel();
    private static final Level testLevel;
    private static final Level ffxLevel;
    static Properties properties;

    static {
        Level lev;
        try {
            lev = Level.parse(System.getProperty("ffx.test.log", "INFO").toUpperCase());
        } catch (Exception ex) {
            logger.warning(String.format(" Exception %s in parsing value of ffx.test.log", ex));
            lev = origLevel;
        }
        testLevel = lev;

        try {
            lev = Level.parse(System.getProperty("ffx.log", "INFO").toUpperCase());
        } catch (Exception ex) {
            logger.warning(String.format(" Exception %s in parsing value of ffx.log", ex));
            lev = origLevel;
        }
        ffxLevel = lev;
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

    /**
     * <p>afterTest.</p>
     */
    @After
    public void afterTest() {
        // All properties are set to the values they were at the beginning of the test.
        System.setProperties(properties);
    }

    /**
     * <p>afterClass.</p>
     */
    @AfterClass
    public static void afterClass() {
        Logger.getLogger("ffx").setLevel(origLevel);
        logger.setLevel(origLevel);
    }
}
