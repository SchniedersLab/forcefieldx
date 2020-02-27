//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
import java.io.IOException;
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
public class HistogramReader extends BufferedReader {

    private static final Logger logger = Logger.getLogger(HistogramReader.class.getName());
    private double temperature;
    private double thetaMass;
    private double thetaFriction;
    private double biasMag;
    private int biasCutoff;
    private int countInterval;
    private int lambdaBins;
    private int FLambdaBins;
    private double minFLambda;
    private double dFL;
    private int temperingFlag;
    private double[][] counts;

    /**
     * Private reference to the Histogram instance, if any.
     */
    private final Histogram histogram;

    public HistogramReader(Reader reader) {
        this(null, reader);
    }
    
    /**
     * Constructor.
     *
     * @param histogram The Histogram instance.
     * @param reader    The Reader to use.
     */
    public HistogramReader(Histogram histogram, Reader reader) {
        super(reader);
        this.histogram = histogram;
    }

    /**
     * Read the histogram file.
     */
    public void readHistogramFile() {
        try {
            temperature = Double.parseDouble(readLine().split(" +")[1]);
            thetaMass = Double.parseDouble(readLine().split(" +")[1]);
            thetaFriction = Double.parseDouble(readLine().split(" +")[1]);
            biasMag = Double.parseDouble(readLine().split(" +")[1]);
            biasCutoff = Integer.parseInt(readLine().split(" +")[1]);
            countInterval = Integer.parseInt(readLine().split(" +")[1]);
            lambdaBins = Integer.parseInt(readLine().split(" +")[1]);
            FLambdaBins = Integer.parseInt(readLine().split(" +")[1]);
            minFLambda = Double.parseDouble(readLine().split(" +")[1]);
            dFL = Double.parseDouble(readLine().split(" +")[1]);
            temperingFlag = Integer.parseInt(readLine().split(" +")[1]);

            counts = new double[lambdaBins][FLambdaBins];
            for (int i = 0; i < lambdaBins; i++) {
                String[] countToks = readLine().split(" +");
                for (int j = 0; j < FLambdaBins; j++) {
                    counts[i][j] = Double.parseDouble(countToks[j]);
                }
            }

            if (histogram != null) {
                applyToHistogram();
            }
        } catch (Exception e) {
            String message = " Invalid OST Histogram file.";
            logger.log(Level.SEVERE, message, e);
        }
        try {
            close();
        } catch (IOException ioe) {
            String histoName = histogram == null ? "unknown file" : histogram.toString();
            logger.warning(String.format(" Failed to close histogram reader for %s", histoName));
        }
    }

    public int getLambdaBins() {
        return lambdaBins;
    }

    public int getFLambdaBins() {
        return FLambdaBins;
    }
    
    private void applyToHistogram() {
        histogram.temperature = temperature;
        histogram.thetaMass = thetaMass;
        histogram.thetaFriction = thetaFriction;
        histogram.biasMag = biasMag;
        histogram.biasCutoff = biasCutoff;
        histogram.countInterval = countInterval;
        histogram.lambdaBins = lambdaBins;
        histogram.FLambda = new double[histogram.lambdaBins];
        histogram.dL = 1.0 / (histogram.lambdaBins - 1);
        histogram.dL_2 = histogram.dL / 2.0;
        histogram.FLambdaBins = FLambdaBins;
        histogram.minFLambda = minFLambda;
        histogram.dFL = dFL;
        histogram.dFL_2 = histogram.dFL / 2.0;
        histogram.setTempering(temperingFlag != 0);

        // Allocate memory for the recursion kernel.
        histogram.allocateRecursionKernel();

        for (int i = 0; i < histogram.lambdaBins; i++) {
            for (int j = 0; j < histogram.FLambdaBins; j++) {
                histogram.setRecursionKernelValue(i, j, counts[i][j]);
            }
        }
    }
}
