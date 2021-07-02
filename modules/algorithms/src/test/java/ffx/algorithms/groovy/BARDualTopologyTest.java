package ffx.algorithms.groovy;

import ffx.algorithms.misc.AlgorithmsTest;
import org.junit.Assert;
import org.junit.Test;

public class BARDualTopologyTest extends AlgorithmsTest {


    /** Tests BAR script with nw input. */
    @Test
    public void testBAR() {
        if (!ffxOpenMM) {
            return;
        }
        
        // Set-up the input arguments for the script.
        String[] args = {"-t", "298","--nw","8","--ac","28", "src/main/java/ffx/algorithms/structures/testBar/dualTopology/HIS_solv.20.xyz", "src/main/java/ffx/algorithms/structures/testBar/dualTopology/HID_solv.20.xyz"};
        binding.setVariable("args", args);
        // Evaluate the script.
        BAR bar = new BAR(binding).run();
        algorithmsScript = bar;

        System.out.println("after run");
        double expectedFepFor= -36.4647;
        double actualFepFor = bar.getFepFor();
        Assert.assertEquals(expectedFepFor, actualFepFor, 0.7);

        double expectedFepBack= -36.0005;
        double actualFepBack = bar.getFepBack();
        Assert.assertEquals(expectedFepBack, actualFepBack, 0.7);

        double expectedhFor= -8.7772;
        double actualhFor = bar.gethFor();
        Assert.assertEquals(expectedhFor, actualhFor, 70);

        double expectedhBack= 27.7404;
        double actualhBack = bar.gethBack();
        Assert.assertEquals(expectedhBack, actualhBack, 60);

        double expectedhBAR= 14.188;
        double actualhBAR = bar.gethBAR();
        Assert.assertEquals(expectedhBAR, actualhBAR, 7);

        double expectedsFor= 0.0929;
        double actualsFor = bar.getsFor();
        Assert.assertEquals(expectedsFor, actualsFor, 0.1);

        double expectedsBack= 0.2139;
        double actualsBack = bar.getsBack();
        Assert.assertEquals(expectedsBack, actualsBack, 0.3);

        double expectedsBAR= 0.1684;
        double actualsBAR = bar.getsBAR();
        Assert.assertEquals(expectedsBAR, actualsBAR, 0.5);

        double expectedBARIteration = -35.7346;
        double actualBARIteration = bar.getBarEnergy();
        Assert.assertEquals(expectedBARIteration, actualBARIteration, 1.0);

        double expectedBARBootstrap= -35.9829;
        double actualBARBootstrap = bar.getBarEnergyBS();
        Assert.assertEquals(expectedBARBootstrap, actualBARBootstrap, 1.0);
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
