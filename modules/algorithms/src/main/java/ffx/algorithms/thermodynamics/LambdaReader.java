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
import static java.lang.String.format;

/**
 * Read in the current value of Lambda, its velocity and the number of
 * counts.
 */
class LambdaReader extends BufferedReader {

    private static final Logger logger = Logger.getLogger(LambdaReader.class.getName());

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
    LambdaReader(TransitionTemperedOSRW transitionTemperedOSRW, Reader reader) {
        super(reader);
        this.transitionTemperedOSRW = transitionTemperedOSRW;
    }

    /**
     * Read the Lambda restart file.
     *
     * @param resetEnergyCount Flag to indicate if the energy count should be read in.
     */
    void readLambdaFile(boolean resetEnergyCount) {
        try {
            transitionTemperedOSRW.lambda = Double.parseDouble(readLine().split(" +")[1]);
            transitionTemperedOSRW.halfThetaVelocity = Double.parseDouble(readLine().split(" +")[1]);
            transitionTemperedOSRW.setLambda(transitionTemperedOSRW.lambda);
        } catch (Exception e) {
            String message = " Invalid OSRW Lambda file.";
            logger.log(Level.SEVERE, message, e);
        }
        if (!resetEnergyCount) {
            try {
                transitionTemperedOSRW.energyCount = Integer.parseUnsignedInt(readLine().split(" +")[1]);
            } catch (Exception e) {
                String message = format(" Could not find number of steps taken in OSRW Lambda file: %s", e.toString());
                logger.log(Level.WARNING, message);
            }
        }
    }
}
