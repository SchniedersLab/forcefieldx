package ffx.utilities;

import org.junit.AfterClass;
import org.junit.BeforeClass;

import java.util.logging.Level;
import java.util.logging.Logger;

public abstract class BaseFFXTest {
    protected static final Logger logger = Logger.getLogger(BaseFFXTest.class.getName());
    private static final Level origLevel = Logger.getLogger("ffx").getLevel();
    private static final Level testLevel;
    private static final Level ffxLevel;

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

    @BeforeClass
    public static void beforeClass() {
        // Set appropriate logging levels for interior/exterior Loggers.
        Logger.getLogger("ffx").setLevel(ffxLevel);
        logger.setLevel(testLevel);
    }

    @AfterClass
    public static void afterClass() {
        Logger.getLogger("ffx").setLevel(origLevel);
        logger.setLevel(origLevel);
    }
}
