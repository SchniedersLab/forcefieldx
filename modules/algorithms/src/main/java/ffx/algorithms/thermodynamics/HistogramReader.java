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
package ffx.algorithms.thermodynamics;

import java.io.BufferedReader;
import java.io.Reader;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Read in the TT-OSRW Histogram.
 */
class HistogramReader extends BufferedReader {

    private static final Logger logger = Logger.getLogger(HistogramReader.class.getName());

    /**
     * Private reference to the TTOSRW instance.
     */
    private TransitionTemperedOSRW transitionTemperedOSRW;

    /**
     * Constructor.
     *
     * @param transitionTemperedOSRW The parent TTOSRW instance.
     * @param reader                 The Reader to use.
     */
    HistogramReader(TransitionTemperedOSRW transitionTemperedOSRW, Reader reader) {
        super(reader);
        this.transitionTemperedOSRW = transitionTemperedOSRW;
    }

    /**
     * Read the histogram file.
     */
    void readHistogramFile() {
        try {
            transitionTemperedOSRW.temperature = Double.parseDouble(readLine().split(" +")[1]);
            transitionTemperedOSRW.thetaMass = Double.parseDouble(readLine().split(" +")[1]);
            transitionTemperedOSRW.thetaFriction = Double.parseDouble(readLine().split(" +")[1]);
            transitionTemperedOSRW.biasMag = Double.parseDouble(readLine().split(" +")[1]);
            transitionTemperedOSRW.biasCutoff = Integer.parseInt(readLine().split(" +")[1]);
            transitionTemperedOSRW.countInterval = Integer.parseInt(readLine().split(" +")[1]);

            transitionTemperedOSRW.lambdaBins = Integer.parseInt(readLine().split(" +")[1]);
            transitionTemperedOSRW.FLambda = new double[transitionTemperedOSRW.lambdaBins];
            transitionTemperedOSRW.dL = 1.0 / (transitionTemperedOSRW.lambdaBins - 1);
            transitionTemperedOSRW.dL_2 = transitionTemperedOSRW.dL / 2.0;

            transitionTemperedOSRW.FLambdaBins = Integer.parseInt(readLine().split(" +")[1]);
            transitionTemperedOSRW.minFLambda = Double.parseDouble(readLine().split(" +")[1]);
            transitionTemperedOSRW.dFL = Double.parseDouble(readLine().split(" +")[1]);
            transitionTemperedOSRW.dFL_2 = transitionTemperedOSRW.dFL / 2.0;

            int flag = Integer.parseInt(readLine().split(" +")[1]);
            transitionTemperedOSRW.setTempering(flag != 0);

            // Allocate memory for the recursion kernel.
            transitionTemperedOSRW.allocateRecursionKernel();

            for (int i = 0; i < transitionTemperedOSRW.lambdaBins; i++) {
                String[] counts = readLine().split(" +");
                for (int j = 0; j < transitionTemperedOSRW.FLambdaBins; j++) {
                    transitionTemperedOSRW.setRecursionKernelValue(i, j, Double.parseDouble(counts[j]));
                }
            }
        } catch (Exception e) {
            String message = " Invalid OSRW Histogram file.";
            logger.log(Level.SEVERE, message, e);
        }
    }
}
