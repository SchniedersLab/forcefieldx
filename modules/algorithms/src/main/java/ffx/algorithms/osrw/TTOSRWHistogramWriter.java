/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.algorithms.osrw;

import java.io.PrintWriter;
import java.io.Writer;

/**
 * Write out the TT-OSRW Histogram.
 */
class TTOSRWHistogramWriter extends PrintWriter {

    /**
     * Private reference to the TTOSRW instance.
     */
    private TransitionTemperedOSRW transitionTemperedOSRW;

    /**
     * Constructor.
     *
     * @param transitionTemperedOSRW The parent TTOSRW instance.
     * @param writer                 The Writer to use.
     */
    TTOSRWHistogramWriter(TransitionTemperedOSRW transitionTemperedOSRW, Writer writer) {
        super(writer);
        this.transitionTemperedOSRW = transitionTemperedOSRW;
    }

    /**
     * Write the Histogram file.
     */
    void writeHistogramFile() {
        printf("Temperature     %15.3f\n", transitionTemperedOSRW.temperature);
        printf("Lambda-Mass     %15.8e\n", transitionTemperedOSRW.thetaMass);
        printf("Lambda-Friction %15.8e\n", transitionTemperedOSRW.thetaFriction);
        printf("Bias-Mag        %15.8e\n", transitionTemperedOSRW.biasMag);
        printf("Bias-Cutoff     %15d\n", transitionTemperedOSRW.biasCutoff);
        printf("Count-Interval  %15d\n", transitionTemperedOSRW.countInterval);
        printf("Lambda-Bins     %15d\n", transitionTemperedOSRW.lambdaBins);
        printf("FLambda-Bins    %15d\n", transitionTemperedOSRW.FLambdaBins);
        printf("Flambda-Min     %15.8e\n", transitionTemperedOSRW.minFLambda);
        printf("Flambda-Width   %15.8e\n", transitionTemperedOSRW.dFL);
        int flag = 0;
        if (transitionTemperedOSRW.isTempering()) {
            flag = 1;
        }
        printf("Tempering       %15d\n", flag);
        for (int i = 0; i < transitionTemperedOSRW.lambdaBins; i++) {
            printf("%g", transitionTemperedOSRW.getRecursionKernelValue(i, 0));
            for (int j = 1; j < transitionTemperedOSRW.FLambdaBins; j++) {
                printf(" %g", transitionTemperedOSRW.getRecursionKernelValue(i, j));
            }
            println();
        }
    }
}
