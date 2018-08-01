package ffx.algorithms;

import edu.rit.pj.Comm;
import ffx.utilities.BaseFFXTest;
import org.junit.BeforeClass;

import java.util.logging.Level;

public abstract class PJDependentTest extends BaseFFXTest {

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
