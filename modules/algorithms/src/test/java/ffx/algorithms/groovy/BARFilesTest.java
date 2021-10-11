//******************************************************************************
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
//******************************************************************************
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
        String[] args = {"-t", "298","--nw","8","--ac","1-13", "--utb", "src/main/java/ffx/algorithms/structures/testBar/singleTopology/dimethylphosphate.100.xyz"};
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
        Assert.assertEquals(expectedhFor, actualhFor, 50);

        double expectedhBack= 123.5150;
        double actualhBack = bar.gethBack();
        Assert.assertEquals(expectedhBack, actualhBack, 90);

        double expectedhBAR= -29.0062;
        double actualhBAR = bar.gethBAR();
        Assert.assertEquals(expectedhBAR, actualhBAR, 4);

        double expectedsFor= -0.1401;
        double actualsFor = bar.getsFor();
        Assert.assertEquals(expectedsFor, actualsFor, 0.15);

        double expectedsBack= 0.8127;
        double actualsBack = bar.getsBack();
        Assert.assertEquals(expectedsBack, actualsBack, 0.3);

        double expectedsBAR= 0.2958;
        double actualsBAR = bar.getsBAR();
        Assert.assertEquals(expectedsBAR, actualsBAR, 0.1);

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
