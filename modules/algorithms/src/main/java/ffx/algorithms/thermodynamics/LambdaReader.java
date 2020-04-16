// ******************************************************************************
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
// ******************************************************************************
package ffx.algorithms.thermodynamics;

import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering.Histogram;
import java.io.BufferedReader;
import java.io.Reader;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Read in the current value of Lambda, its velocity and the number of counts.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class LambdaReader extends BufferedReader {

  private static final Logger logger = Logger.getLogger(LambdaReader.class.getName());
  private double lambda;
  private double halfThetaVel;
  private int nSteps = 0;
  private int histoIndex = 0;
  private boolean resetEnergyCount = false;

  /**
   * Constructor.
   *
   * @param reader The Reader to use.
   */
  public LambdaReader(Reader reader) {
    super(reader);
  }

  /**
   * Get the index of the histogram associated with this lambda file.
   *
   * @return Associated histogram index.
   */
  public int getHistogramIndex() {
    return histoIndex;
  }

  /**
   * Read the Lambda restart file.
   *
   * @param resetEnergyCount Flag to indicate if the energy count should be read in.
   */
  public void readLambdaFile(boolean resetEnergyCount) {
    this.resetEnergyCount = resetEnergyCount;
    try {
      lambda = Double.parseDouble(readLine().split(" +")[1]);
      halfThetaVel = Double.parseDouble(readLine().split(" +")[1]);
      for (int i = 0; i < 2; i++) {
        String line = readLine();
        if (line == null) {
          break;
        }
        if (line.startsWith("Steps-Taken")) {
          nSteps = Integer.parseInt(line.split(" +")[1]);
        } else if (line.startsWith("Histogram")) {
          histoIndex = Integer.parseInt(line.split(" +")[1]);
        }
      }
    } catch (Exception e) {
      String message = " Invalid OST Lambda file.";
      logger.log(Level.SEVERE, message, e);
    }
  }

  void setVariables(OrthogonalSpaceTempering ost) {
    ost.setLambda(lambda);
    Histogram histo = ost.getHistogram();
    histo.setHalfThetaVelocity(halfThetaVel);
    if (!resetEnergyCount && nSteps > 0) {
      ost.setEnergyCount(nSteps);
    }
  }
}
