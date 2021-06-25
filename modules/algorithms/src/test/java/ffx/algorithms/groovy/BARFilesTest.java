package ffx.algorithms.groovy;

import ffx.algorithms.misc.AlgorithmsTest;
import org.junit.Assert;
import org.junit.Test;

public class BARFilesTest extends AlgorithmsTest {


    /** Tests BAR script with tinker bar files input. */
    @Test
    public void testBARFiles() {

        if (!ffxOpenMM) {
            return;
        }
        
        // Set-up the input arguments for the script.
        String[] args = {"-t", "298","--nw","8","--ac","1-13", "--utb", "src/main/java/ffx/algorithms/structures/testBar/dimethylphosphate.100.xyz"};
        binding.setVariable("args", args);
        // Evaluate the script.
        BAR bar = new BAR(binding).run();
        algorithmsScript = bar;

        double expectedFepFor= -115.4001;
        double actualFepFor = bar.getFepFor();
        Assert.assertEquals(expectedFepFor, actualFepFor, 0.5);

        double expectedFepBack= -118.6607;
        double actualFepBack = bar.getFepBack();
        Assert.assertEquals(expectedFepBack, actualFepBack, 0.5);

        double expectedhFor= -157.1466;
        double actualhFor = bar.gethFor();
        Assert.assertEquals(expectedhFor, actualhFor, 33);

        double expectedhBack= 123.5150;
        double actualhBack = bar.gethBack();
        Assert.assertEquals(expectedhBack, actualhBack, 70);

        double expectedhBAR= -29.0062;
        double actualhBAR = bar.gethBAR();
        Assert.assertEquals(expectedhBAR, actualhBAR, 6);

        double expectedsFor= -0.1401;
        double actualsFor = bar.getsFor();
        Assert.assertEquals(expectedsFor, actualsFor, 0.05);

        double expectedsBack= 0.8127;
        double actualsBack = bar.getsBack();
        Assert.assertEquals(expectedsBack, actualsBack, 0.05);

        double expectedsBAR= 0.2958;
        double actualsBAR = bar.getsBAR();
        Assert.assertEquals(expectedsBAR, actualsBAR, 0.05);

        double expectedBARIteration = -117.1754;
        double actualBARIteration = bar.getBarEnergy();
        Assert.assertEquals(expectedBARIteration, actualBARIteration, 0.5);

        double expectedBARBootstrap= -117.1464;
        double actualBARBootstrap = bar.getBarEnergyBS();
        Assert.assertEquals(expectedBARBootstrap, actualBARBootstrap, 0.3);
    }
    @Test
    public void testBARHelp() {
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        BAR bar = new BAR(binding).run();
        algorithmsScript = bar;
    }

}
