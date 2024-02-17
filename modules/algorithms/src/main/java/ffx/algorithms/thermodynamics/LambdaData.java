// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import ffx.utilities.FileUtils;

import javax.annotation.Nullable;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessOrder;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorOrder;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import java.io.File;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.asin;
import static org.apache.commons.math3.util.FastMath.sqrt;

@XmlRootElement(name = "LambdaData")
@XmlAccessorOrder(XmlAccessOrder.ALPHABETICAL)
@XmlAccessorType(XmlAccessType.NONE)
public class LambdaData {

  private static final Logger logger = Logger.getLogger(LambdaData.class.getName());

  /**
   * State variable lambda ranges from 0.0 .. 1.0.
   */
  @XmlElement(name = "Lambda", defaultValue = "0.0")
  double lambda = 0.0;

  /**
   * Velocity of the theta particle.
   */
  @XmlElement(name = "ThetaVelocity", defaultValue = "0.0")
  double thetaVelocity = 0.0;

  /**
   * Velocity of the theta particle.
   */
  @XmlElement(name = "ThetaAcceleration", defaultValue = "0.0")
  double thetaAcceleration = 0.0;

  /**
   * The number of lambda steps taken by the propagate lambda method.
   */
  @XmlElement(name = "StepsTaken", defaultValue = "0")
  long stepsTaken = 0;

  /**
   * Index of the current Histogram.
   */
  @XmlElement(name = "HistogramIndex", defaultValue = "0")
  int histogramIndex = 0;

  /**
   * The LambdaData file to write to / read from.
   */
  private File lambdaFile = null;

  /**
   * The LambdaData file name.
   */
  private String lambdaFileName = null;

  /**
   * Flag indicating if this LambdaData came from reading a file.
   */
  private boolean lambdaRead = false;

  /**
   * Map lambda to a periodic variable theta.
   * <code>theta = asin(sqrt(lambda))</code>
   * <code>lambda = sin^2 (theta).</code>
   */
  double theta = 0.0;

  private void rectify() {
    theta = asin(sqrt(lambda));
  }

  /**
   * Marshall the lambda data to a file.
   *
   * @param lambdaData The LambdaData to write.
   * @param lambdaFile The file to write to.
   */
  public static void writeLambdaData(LambdaData lambdaData, File lambdaFile) {
    try {
      JAXBContext jaxbContext = JAXBContext.newInstance(LambdaData.class);
      Marshaller jaxbMarshaller = jaxbContext.createMarshaller();
      jaxbMarshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
      jaxbMarshaller.marshal(lambdaData, lambdaFile);
    } catch (Exception e) {
      logger.warning(" Exception writing histogram:\n " + e);
    }
  }

  /**
   * Unmarshall the histogram data.
   *
   * @param file The file to read from.
   * @return The Histogram data.
   */
  public static LambdaData readLambdaData(@Nullable File file) {
    if (file != null && file.exists() && file.canRead()) {
      try {
        JAXBContext jaxbContext = JAXBContext.newInstance(LambdaData.class);
        Unmarshaller unmarshaller = jaxbContext.createUnmarshaller();
        LambdaData lambdaData = (LambdaData) unmarshaller.unmarshal(file);
        lambdaData.setLambdaFile(file);
        lambdaData.setLambdaRead(true);
        lambdaData.rectify();
        return lambdaData;
      } catch (Exception e) {
        logger.warning(" Exception reading LambdaData:\n " + e);
      }
    }
    LambdaData lambdaData = new LambdaData();
    lambdaData.setLambdaFile(file);
    lambdaData.rectify();
    return lambdaData;
  }

  /**
   * Marshall the lambda data to a file.
   */
  public void writeLambdaData() {
    writeLambdaData(this, lambdaFile);
  }

  public File getLambdaFile() {
    return lambdaFile;
  }

  public void setLambdaFile(@Nullable File lambdaFile) {
    this.lambdaFile = lambdaFile;
    if (lambdaFile != null) {
      lambdaFileName = FileUtils.relativePathTo(lambdaFile).toString();
    } else {
      lambdaFileName = null;
    }
  }

  public String getLambdaFileName() {
    return lambdaFileName;
  }

  public boolean wasLambdaRead() {
    return lambdaRead;
  }

  public void setLambdaRead(boolean lambdaRead) {
    this.lambdaRead = lambdaRead;
  }

  public void setHistogramIndex(int histogramIndex) {
    this.histogramIndex = histogramIndex;
  }

  public void setStepsTaken(int stepsTaken) {
    this.stepsTaken = stepsTaken;
  }
}
