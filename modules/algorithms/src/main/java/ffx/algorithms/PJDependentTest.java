package ffx.algorithms;

import java.util.logging.Level;

import org.junit.BeforeClass;

import edu.rit.pj.Comm;

import ffx.utilities.BaseFFXTest;

/**
 * <p>Abstract PJDependentTest class.</p>
 *
 * @author Michael J. Schnieders
 */
public abstract class PJDependentTest extends BaseFFXTest {

    /**
     * <p>beforeClass.</p>
     */
    @BeforeClass
    public static void beforeClass() {
        BaseFFXTest.beforeClass();
        // Initialize Parallel Java
        try {
            Comm.world();
        } catch (IllegalStateException ise) {
            try {
                String args[] = new String[0];
                Comm.init(args);
            } catch (Exception e) {
                String message = " Exception starting up the Parallel Java communication layer.";
                logger.log(Level.WARNING, message, e.toString());
                message = " Skipping rotamer optimization test.";
                logger.log(Level.WARNING, message, e.toString());
            }
        }
    }
}
