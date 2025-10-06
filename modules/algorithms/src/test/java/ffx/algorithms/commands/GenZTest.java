package ffx.algorithms.commands;

import ffx.algorithms.misc.AlgorithmsTest;
import org.junit.Assert;
import org.junit.Test;

public class GenZTest extends AlgorithmsTest{

    @Test
    public void testGenZ(){
        if (!ffxOpenMM) {
            return;
        }

        String filepath = getResourcePath("DEHK.genz.pdb");
        String[] args = {"--tR", "--pH", "7.0", "--rEE", "2.0", "--kPH 2.0", "-O","--pKa", filepath};
        binding.setVariable("args", args);

        GenZ genZ = new GenZ(binding).run();
        algorithmsScript = genZ;

        double expectedAspPop = 1.00;
        double actualAspPop = genZ.getPopulationArray()[0][0];
        Assert.assertEquals(expectedAspPop, actualAspPop, 0.0001);

        double expectedAshPop = 0.00;
        double actualAshPop = genZ.getPopulationArray()[0][1];
        Assert.assertEquals(expectedAshPop, actualAshPop, 0.0001);

        double expectedGluPop = 1.00;
        double actualGluPop = genZ.getPopulationArray()[1][0];
        Assert.assertEquals(expectedGluPop, actualGluPop, 0.0001);

        double expectedGlhPop = 0.00;
        double actualGlhPop = genZ.getPopulationArray()[1][1];
        Assert.assertEquals(expectedGlhPop, actualGlhPop, 0.0001);

        double expectedHisPop = 0.00;
        double actualHisPop = genZ.getPopulationArray()[2][0];
        Assert.assertEquals(expectedHisPop, actualHisPop, 0.001);

        double expectedHiePop = 0.285801;
        double actualHiePop = genZ.getPopulationArray()[2][1];
        Assert.assertEquals(expectedHiePop, actualHiePop, 0.001);

        double expectedHidPop = 0.714199;
        double actualHidPop = genZ.getPopulationArray()[2][1];
        Assert.assertEquals(expectedHidPop, actualHidPop, 0.001);

        double expectedLysPop = 1.00;
        double actualLysPop = genZ.getPopulationArray()[3][0];
        Assert.assertEquals(expectedLysPop, actualLysPop, 0.0001);

        double expectedLydPop = 0.00;
        double actualLydPop = genZ.getPopulationArray()[3][1];
        Assert.assertEquals(expectedLydPop, actualLydPop, 0.0001);

    }

    @Test
    public void testGenZHelp() {
        String[] args = {"-h"};
        binding.setVariable("args", args);

        // Evaluate the script.
        GenZ genZ = new GenZ(binding).run();
        algorithmsScript = genZ;
    }
}
