package ffx.algorithms.groovy;

import ffx.algorithms.misc.AlgorithmsTest;
import org.junit.Assert;
import org.junit.Test;

public class BARTest extends AlgorithmsTest {


    /** Tests BAR script with nw input. */
    @Test
    public void testBAR() {
        if (!ffxOpenMM) {
            return;
        }
        
        // Set-up the input arguments for the script.
        String[] args = {"--t1", "298","--nw","8","--ac","1-13", "src/main/java/ffx/algorithms/structures/testBar/dimethylphosphate.10.xyz"};
        binding.setVariable("args", args);
        // Evaluate the script.
        BAR bar = new BAR(binding).run();
        algorithmsScript = bar;

        System.out.println("after run");
        double expectedFepFor= -115.2142;
        double actualFepFor = bar.getFepFor();
        Assert.assertEquals(expectedFepFor, actualFepFor, 0.2);

        double expectedFepBack= -118.4924;
        double actualFepBack = bar.getFepBack();
        Assert.assertEquals(expectedFepBack, actualFepBack, 0.2);

        double expectedhFor= -32.0327;
        double actualhFor = bar.gethFor();
        Assert.assertEquals(expectedhFor, actualhFor, 7);

        double expectedhBack= 52.9186;
        double actualhBack = bar.gethBack();
        Assert.assertEquals(expectedhBack, actualhBack, 0.1);

        double expectedhBAR= 14.6302;
        double actualhBAR = bar.gethBAR();
        Assert.assertEquals(expectedhBAR, actualhBAR, 1);

        double expectedsFor= 0.0771;
        double actualsFor = bar.getsFor();
        Assert.assertEquals(expectedsFor, actualsFor, 0.05);

        double expectedsBack= 0.3552;
        double actualsBack = bar.getsBack();
        Assert.assertEquals(expectedsBack, actualsBack, 0.001);

        double expectedsBAR= 0.2280;
        double actualsBAR = bar.getsBAR();
        Assert.assertEquals(expectedsBAR, actualsBAR, 0.01);

        double expectedBARIteration = -53.3234;
        double actualBARIteration = bar.getBarEnergy();
        Assert.assertEquals(expectedBARIteration, actualBARIteration, 0.1);

        double expectedBARBootstrap= -53.3224;
        double actualBARBootstrap = bar.getBarEnergyBS();
        Assert.assertEquals(expectedBARBootstrap, actualBARBootstrap, 0.2);
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
