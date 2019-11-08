//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.algorithms.thermodynamics;

import java.io.BufferedReader;
import java.io.Reader;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering.Histogram;

/**
 * Read in the Histogram.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
class HistogramReader extends BufferedReader {

    private static final Logger logger = Logger.getLogger(HistogramReader.class.getName());

    /**
     * Private reference to the Histogram instance.
     */
    private Histogram histogram;

    /**
     * Constructor.
     *
     * @param histogram The Histogram instance.
     * @param reader    The Reader to use.
     */
    HistogramReader(Histogram histogram, Reader reader) {
        super(reader);
        this.histogram = histogram;
    }

    /**
     * Read the histogram file.
     */
    void readHistogramFile() {
        try {
            histogram.temperature = Double.parseDouble(readLine().split(" +")[1]);
            histogram.thetaMass = Double.parseDouble(readLine().split(" +")[1]);
            histogram.thetaFriction = Double.parseDouble(readLine().split(" +")[1]);
            histogram.biasMag = Double.parseDouble(readLine().split(" +")[1]);
            histogram.biasCutoff = Integer.parseInt(readLine().split(" +")[1]);
            histogram.countInterval = Integer.parseInt(readLine().split(" +")[1]);
            histogram.lambdaBins = Integer.parseInt(readLine().split(" +")[1]);
            histogram.FLambda = new double[histogram.lambdaBins];
            histogram.dL = 1.0 / (histogram.lambdaBins - 1);
            histogram.dL_2 = histogram.dL / 2.0;
            histogram.FLambdaBins = Integer.parseInt(readLine().split(" +")[1]);
            histogram.minFLambda = Double.parseDouble(readLine().split(" +")[1]);
            histogram.dFL = Double.parseDouble(readLine().split(" +")[1]);
            histogram.dFL_2 = histogram.dFL / 2.0;

            int flag = Integer.parseInt(readLine().split(" +")[1]);
            histogram.setTempering(flag != 0);

            // Allocate memory for the recursion kernel.
            histogram.allocateRecursionKernel();

            for (int i = 0; i < histogram.lambdaBins; i++) {
                String[] counts = readLine().split(" +");
                for (int j = 0; j < histogram.FLambdaBins; j++) {
                    histogram.setRecursionKernelValue(i, j, Double.parseDouble(counts[j]));
                }
            }
        } catch (Exception e) {
            String message = " Invalid OST Histogram file.";
            logger.log(Level.SEVERE, message, e);
        }
    }
}
