//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.numerics.estimator;
import ffx.utilities.Constants;
import ffx.utilities.FFXTest;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

public class MBARHarmonicOscillatorsTest extends FFXTest {

    /**
     * Test the MBAR estimator numerics with harmonic oscillators. This test uses SCI & NR, so it
     * isn't exactly deterministic. The error is set to 1.0E-1, so the first decimal is fine for this.
     * Pymbar does the same thing.
     */
    @Test
    public void testMBAROscillators() {
        double[] O_k = {1, 2, 3, 4};
        double[] K_k = {.5, 1.0, 1.5, 2};
        int[] N_k = {10000, 10000, 10000, 10000}; // No support for different number of snapshots
        double beta = 1.0;

        // Create an instance of HarmonicOscillatorsTestCase
        MultistateBennettAcceptanceRatio.HarmonicOscillatorsTestCase testCase = new MultistateBennettAcceptanceRatio.HarmonicOscillatorsTestCase(O_k, K_k, beta);

        // Generate sample data
        String setting = "u_kln";
        Object[] sampleResult = testCase.sample(N_k, setting, (long) 0);
        double[][][] u_kln = (double[][][]) sampleResult[1];
        double[] temps = {1 / Constants.R};

        MultistateBennettAcceptanceRatio mbar = new MultistateBennettAcceptanceRatio(O_k, u_kln, temps, 1.0E-7, MultistateBennettAcceptanceRatio.SeedType.ZEROS);
        double[] mbarFEEstimates = mbar.getMBARFreeEnergies();
        double[] mbarErrorEstimates = mbar.getBinUncertainties();
        double[][] mbarDiffMatrix = mbar.getDiffMatrix();
        double[] mbarFEExpected = new double[]{0.0, 0.3468272332334239, 0.554882810046907, 0.6909139007747198};
        double[] mbarErrorExpected = new double[]{0.00647778279366289, 0.006176323555016366, 0.008170508071621832};
        double[][] mbarDiffMatrixExpected = new double[][]{
                {0.0, 0.00647778279366289, 0.010874539771152386, 0.01568591641036036},
                {0.005859375, 0.0, 0.006176323555016366, 0.012881744099875898},
                {0.010697706201272776, 0.00647778279366289, 0.0, 0.008170508071621832},
                {0.015563845166512918, 0.013028968812623373, 0.008170508071621832, 0.0}
        };

        // Compare tolerance to pymbar --> I tested this with a for loop over 1000 times and it always passed
        for (int i = 0; i < mbarFEExpected.length; i++) {
            Assert.assertEquals(mbarFEExpected[i], mbarFEEstimates[i], 1.0E-1);
        }
        for (int i = 0; i < mbarErrorExpected.length; i++) {
            Assert.assertEquals(mbarErrorExpected[i], mbarErrorEstimates[i], 1.0E-1);
        }
        for (int i = 0; i < mbarDiffMatrixExpected.length; i++) {
            for (int j = 0; j < mbarDiffMatrixExpected[i].length; j++) {
                Assert.assertEquals(mbarDiffMatrixExpected[i][j], mbarDiffMatrix[i][j], 1.0E-1);
            }
        }
    }
}
