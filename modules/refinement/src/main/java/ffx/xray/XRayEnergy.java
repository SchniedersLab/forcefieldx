// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.xray;

import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parameters.ForceField;
import ffx.xray.refine.RefinementMode;
import ffx.xray.refine.RefinementModel;

import javax.annotation.Nullable;
import java.util.List;
import java.util.logging.Logger;

import static ffx.numerics.math.MatrixMath.mat3Determinant;
import static ffx.numerics.math.ScalarMath.b2u;
import static ffx.numerics.math.ScalarMath.u2b;
import static ffx.utilities.Constants.R;
import static java.lang.Math.pow;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static org.apache.commons.math3.util.FastMath.PI;

/**
 * Combine the X-ray target and chemical potential energy.
 *
 * @author Timothy D. Fenn
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class XRayEnergy implements LambdaInterface, CrystalPotential {

  private static final Logger logger = Logger.getLogger(XRayEnergy.class.getName());
  private static final double eightPI2 = 8.0 * PI * PI;
  private static final double eightPI23 = eightPI2 * eightPI2 * eightPI2;

  private final DiffractionData diffractionData;
  private final RefinementModel refinementModel;
  private final RefinementMode refinementMode;
  private final Atom[] activeAtomArray;
  private double[] optimizationScaling = null;
  private final double kTbNonzero;
  private final double kTbSimWeight;
  private final boolean lambdaTerm;
  private final double[] g2;
  private final double[] dUdXdL;
  protected double lambda = 1.0;
  private final int nXYZ;
  private final int nB;
  private final int nOCC;
  private boolean xrayTerms = true;
  private boolean restraintTerms = true;
  private double totalEnergy;
  private double dEdL;
  private STATE state = STATE.BOTH;

  /**
   * Diffraction data energy target
   *
   * @param diffractionData {@link ffx.xray.DiffractionData} object to associate with the target
   */
  public XRayEnergy(DiffractionData diffractionData) {
    this.diffractionData = diffractionData;

    refinementModel = diffractionData.getRefinementModel();
    refinementMode = refinementModel.getRefinementMode();
    nXYZ = refinementModel.getNumCoordParameters();
    nB = refinementModel.getNumBFactorParameters();
    nOCC = refinementModel.getNumOccupancyParameters();

    double temperature = 50.0;
    kTbNonzero = R * temperature * diffractionData.getbNonZeroWeight();
    kTbSimWeight = R * temperature * diffractionData.getbSimWeight();

    ForceField forceField = diffractionData.getAssembly()[0].getForceField();
    lambdaTerm = forceField.getBoolean("LAMBDATERM", false);

    activeAtomArray = refinementModel.getActiveAtoms();
    int count = activeAtomArray.length;
    dUdXdL = new double[count * 3];
    g2 = new double[count * 3];

    if (refinementMode.includesBFactors()) {
      logger.info("\n B-Factor Refinement Parameters");
      logger.info("  Temperature:                 " + temperature);
      logger.info("  Non-zero restraint weight:   " + diffractionData.getbNonZeroWeight());
      logger.info("  Similarity restraint weight: " + diffractionData.getbSimWeight());
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean destroy() {
    return diffractionData.destroy();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energy(double[] x) {
    double e = 0.0;

    // Unscale the coordinates.
    unscaleCoordinates(x);

    // Load the parameters.
    refinementModel.setParameters(x);

    // Update the coordinates for the experimental term if coordinates are being refined.
    if (refinementMode.includesCoordinates()) {
      diffractionData.updateCoordinates();
    }

    if (xrayTerms) {

      if (lambdaTerm) {
        diffractionData.setLambdaTerm(false);
      }

      // Compute new structure factors.
      diffractionData.computeAtomicDensity();
      // Compute crystal likelihood.
      e = diffractionData.computeLikelihood();

      if (lambdaTerm) {

        // Turn off all atoms scaled by lambda.
        diffractionData.setLambdaTerm(true);

        // Compute new structure factors.
        diffractionData.computeAtomicDensity();

        // Compute crystal likelihood.
        double e2 = diffractionData.computeLikelihood();

        dEdL = e - e2;

        e = lambda * e + (1.0 - lambda) * e2;

        diffractionData.setLambdaTerm(false);
      }
    }

    if (restraintTerms) {
      if (refinementMode.includesBFactors()) {
        // add B restraints
        e += getBFactorRestraints(false);
      }
    }

    // Scale the coordinates and gradients.
    scaleCoordinates(x);

    totalEnergy = e;
    return e;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energyAndGradient(double[] x, double[] g) {
    double e = 0.0;

    // Unscale the coordinates.
    unscaleCoordinates(x);

    // Set the parameters.
    refinementModel.setParameters(x);

    // Zero out the gradient for each atom that holds refinement parameters.
    refinementModel.zeroGradient();

    // Update the coordinates for the experimental term if coordinates are being refined.
    if (refinementMode.includesCoordinates()) {
      diffractionData.updateCoordinates();
    }

    if (xrayTerms) {

      if (lambdaTerm) {
        diffractionData.setLambdaTerm(false);
      }

      // compute new structure factors
      diffractionData.computeAtomicDensity();

      // compute crystal likelihood
      e = diffractionData.computeLikelihood();

      // Compute the X-ray gradients.
      diffractionData.computeAtomicGradients(refinementMode);

      if (lambdaTerm) {
        logger.severe(" Lambda Refinement is not supported.");
        int n = dUdXdL.length;
        arraycopy(g, 0, dUdXdL, 0, n);

        for (Atom a : activeAtomArray) {
          a.setXYZGradient(0.0, 0.0, 0.0);
          a.setLambdaXYZGradient(0.0, 0.0, 0.0);
        }

        // Turn off all atoms scaled by lambda.
        diffractionData.setLambdaTerm(true);

        // Compute new structure factors.
        diffractionData.computeAtomicDensity();

        // Compute crystal likelihood.
        double e2 = diffractionData.computeLikelihood();

        // compute the crystal gradients
        diffractionData.computeAtomicGradients(refinementMode);

        dEdL = e - e2;
        e = lambda * e + (1.0 - lambda) * e2;
        // getXYZGradients(g2);
        for (int i = 0; i < g.length; i++) {
          dUdXdL[i] -= g2[i];
          g[i] = lambda * g[i] + (1.0 - lambda) * g2[i];
        }

        diffractionData.setLambdaTerm(false);
      }
    }

    if (restraintTerms) {
      if (refinementMode.includesBFactors()) {
        // Add B-factor restraints.
        e += getBFactorRestraints(true);
      }
    }

    // Load the gradient over all parameters.
    refinementModel.getGradient(g);

    // Scale the coordinates and gradients.
    scaleCoordinatesAndGradient(x, g);

    totalEnergy = e;
    return e;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getCoordinates(double[] x) {
    if (x == null || x.length != refinementModel.getNumParameters()) {
      x = new double[refinementModel.getNumParameters()];
    }
    refinementModel.getParameters(x);
    return x;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setCoordinates(double[] x) {
    refinementModel.setParameters(x);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public Crystal getCrystal() {
    return diffractionData.getCrystal()[0];
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setCrystal(Crystal crystal) {
    logger.severe(" XRayEnergy does implement setCrystal yet.");
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public STATE getEnergyTermState() {
    return state;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setEnergyTermState(STATE state) {
    this.state = state;
    switch (state) {
      case FAST:
        xrayTerms = false;
        restraintTerms = true;
        break;
      case SLOW:
        xrayTerms = true;
        restraintTerms = false;
        break;
      default:
        xrayTerms = true;
        restraintTerms = true;
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getLambda() {
    return lambda;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setLambda(double lambda) {
    if (lambda <= 1.0 && lambda >= 0.0) {
      this.lambda = lambda;
    } else {
      String message = format("Lambda value %8.3f is not in the range [0..1].", lambda);
      logger.warning(message);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getMass() {
    double[] mass = new double[nXYZ + nB + nOCC];
    refinementModel.getMass(mass);
    return mass;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfVariables() {
    return nXYZ + nB + nOCC;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getScaling() {
    return optimizationScaling;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setScaling(@Nullable double[] scaling) {
    optimizationScaling = scaling;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalEnergy() {
    return totalEnergy;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Return a reference to each variables type.
   */
  @Override
  public VARIABLE_TYPE[] getVariableTypes() {
    VARIABLE_TYPE[] vtypes = new VARIABLE_TYPE[nXYZ + nB + nOCC];
    int i = 0;
    if (refinementMode.includesCoordinates()) {
      for (Atom a : activeAtomArray) {
        vtypes[i++] = VARIABLE_TYPE.X;
        vtypes[i++] = VARIABLE_TYPE.Y;
        vtypes[i++] = VARIABLE_TYPE.Z;
      }
    }
    if (refinementMode.includesBFactors()) {
      for (int j = i; j < nXYZ + nB; i++, j++) {
        vtypes[j] = VARIABLE_TYPE.OTHER;
      }
    }
    if (refinementMode.includesOccupancies()) {
      for (int j = i; j < nXYZ + nB + nOCC; i++, j++) {
        vtypes[j] = VARIABLE_TYPE.OTHER;
      }
    }
    return vtypes;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getVelocity(double[] velocity) {
    if (velocity == null || velocity.length != refinementModel.getNumParameters()) {
      velocity = new double[refinementModel.getNumParameters()];
    }
    refinementModel.getVelocity(velocity);
    return velocity;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setVelocity(double[] velocity) {
    refinementModel.setVelocity(velocity);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getd2EdL2() {
    return 0.0;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getdEdL() {
    return dEdL;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void getdEdXdL(double[] gradient) {
    int n = dUdXdL.length;
    arraycopy(dUdXdL, 0, gradient, 0, n);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setAcceleration(double[] acceleration) {
    refinementModel.setAcceleration(acceleration);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getAcceleration(double[] acceleration) {
    if (acceleration == null || acceleration.length != refinementModel.getNumParameters()) {
      acceleration = new double[refinementModel.getNumParameters()];
    }
    refinementModel.getAcceleration(acceleration);
    return acceleration;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setPreviousAcceleration(double[] previousAcceleration) {
    refinementModel.setPreviousAcceleration(previousAcceleration);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getPreviousAcceleration(double[] previousAcceleration) {
    if (previousAcceleration == null || previousAcceleration.length != refinementModel.getNumParameters()) {
      previousAcceleration = new double[refinementModel.getNumParameters()];
    }
    refinementModel.getPreviousAcceleration(previousAcceleration);
    return previousAcceleration;
  }

  /**
   * determine similarity and non-zero B factor restraints (done independently of
   * getBFactorGradients), affects atomic gradients
   *
   * @param gradient Compute the gradient.
   * @return The energy of the B-factor restraints.
   */
  private double getBFactorRestraints(boolean gradient) {
    double e = 0.0;
    double[] anisou1 = new double[6];
    double[] anisou2;
    double[] gradu = new double[6];
    double threeHalves = 3.0 / 2.0;
    double oneHalf = 1.0 / 2.0;

    // Single B-Factor Restraint.
    for (Atom a : activeAtomArray) {
      if (a.getAnisou(null) == null) {
        // Non-zero restraint using ADP as the partition function.
        // The constant offset below is neglected in the potential.
        // -T S = kB T * ln(Biso^3/2) + ln(C)
        // -T dS = kB T * 3/2 / Biso
        double biso = a.getTempFactor();
        e += -kTbNonzero * Math.log(pow(biso, threeHalves));
        if (gradient) {
          double gradb = -kTbNonzero * threeHalves / biso;
          a.addToTempFactorGradient(gradb);
        }
      } else {
        // Anisotropic B restraint
        anisou1 = a.getAnisou(anisou1);
        // Non-zero restraint using ADP as the partition function.
        // -T S     = -1/2 kB T * ln [(2 PI e)^2 det(U)]
        //          = -1/2 kB T * ln [det(U)] + C
        // -T dS/dU = -1/2 kB T * U^-1
        double det = mat3Determinant(anisou1);
        e += u2b(-oneHalf * kTbNonzero * Math.log(det));
        if (gradient) {
          gradu[0] = u2b(-oneHalf * kTbNonzero * ((anisou1[1] * anisou1[2] - anisou1[5] * anisou1[5]) / det));
          gradu[1] = u2b(-oneHalf * kTbNonzero * ((anisou1[0] * anisou1[2] - anisou1[4] * anisou1[4]) / det));
          gradu[2] = u2b(-oneHalf * kTbNonzero * ((anisou1[0] * anisou1[1] - anisou1[3] * anisou1[3]) / det));
          gradu[3] = u2b(-oneHalf * kTbNonzero * ((2.0 * (anisou1[4] * anisou1[5] - anisou1[3] * anisou1[2])) / det));
          gradu[4] = u2b(-oneHalf * kTbNonzero * ((2.0 * (anisou1[3] * anisou1[5] - anisou1[4] * anisou1[1])) / det));
          gradu[5] = u2b(-oneHalf * kTbNonzero * ((2.0 * (anisou1[3] * anisou1[4] - anisou1[5] * anisou1[0])) / det));
          a.addToAnisouGradient(gradu);
        }
      }
    }

    // B-Factor Restraints.
    List<Atom[]> bonds = refinementModel.getBFactorRestraints();
    for (Atom[] atoms : bonds) {
      Atom a1 = atoms[0];
      Atom a2 = atoms[1];
      boolean isAnisou1 = a1.getAnisou(null) != null;
      boolean isAnisou2 = a2.getAnisou(null) != null;
      if (!isAnisou1 && !isAnisou2) {
        // Both atoms are isotropic.
        double b1 = a1.getTempFactor();
        double b2 = a2.getTempFactor();
        double bdiff = b1 - b2;
        e += kTbSimWeight * bdiff * bdiff;
        if (gradient) {
          double gradb = 2.0 * kTbSimWeight * bdiff;
          a1.addToTempFactorGradient(gradb);
          a2.addToTempFactorGradient(-gradb);
        }
      } else if (isAnisou1 && isAnisou2) {
        // Both atoms are anisotropic.
        anisou1 = a1.getAnisou(anisou1);
        anisou2 = a2.getAnisou(anisou1);
        double det1 = mat3Determinant(anisou1);
        double det2 = mat3Determinant(anisou2);
        double bdiff = det1 - det2;
        double bdiff2 = bdiff * bdiff;
        e += eightPI23 * kTbSimWeight * bdiff2;
        if (gradient) {
          double gradb = eightPI23 * 2.0 * kTbSimWeight * bdiff;
          // Atom 1
          gradu[0] = gradb * (anisou1[1] * anisou1[2] - anisou1[5] * anisou1[5]);
          gradu[1] = gradb * (anisou1[0] * anisou1[2] - anisou1[4] * anisou1[4]);
          gradu[2] = gradb * (anisou1[0] * anisou1[1] - anisou1[3] * anisou1[3]);
          gradu[3] = gradb * (2.0 * (anisou1[4] * anisou1[5] - anisou1[3] * anisou1[2]));
          gradu[4] = gradb * (2.0 * (anisou1[3] * anisou1[5] - anisou1[4] * anisou1[1]));
          gradu[5] = gradb * (2.0 * (anisou1[3] * anisou1[4] - anisou1[5] * anisou1[0]));
          a1.addToAnisouGradient(gradu);
          // Atom 2
          gradu[0] = gradb * (anisou2[5] * anisou2[5] - anisou2[1] * anisou2[2]);
          gradu[1] = gradb * (anisou2[4] * anisou2[4] - anisou2[0] * anisou2[2]);
          gradu[2] = gradb * (anisou2[3] * anisou2[3] - anisou2[0] * anisou2[1]);
          gradu[3] = gradb * (2.0 * (anisou2[3] * anisou2[2] - anisou2[4] * anisou2[5]));
          gradu[4] = gradb * (2.0 * (anisou2[4] * anisou2[1] - anisou2[3] * anisou2[5]));
          gradu[5] = gradb * (2.0 * (anisou2[5] * anisou2[0] - anisou2[3] * anisou2[4]));
          a2.addToAnisouGradient(gradu);
        }
      } else {
        if (!isAnisou1) {
          // Swap a1 and a2.
          a1 = atoms[1];
          a2 = atoms[0];
        }
        anisou1 = a1.getAnisou(anisou1);
        double u2 = b2u(a2.getTempFactor());
        double det1 = mat3Determinant(anisou1);
        // Determinant of a diagonal matrix.
        double det2 = u2 * u2 * u2;
        double bdiff = det1 - det2;
        double bdiff2 = bdiff * bdiff;
        e += eightPI23 * kTbSimWeight * bdiff2;
        if (gradient) {
          double gradb = eightPI23 * 2.0 * kTbSimWeight * bdiff;
          // Atom 1
          gradu[0] = gradb * (anisou1[1] * anisou1[2] - anisou1[5] * anisou1[5]);
          gradu[1] = gradb * (anisou1[0] * anisou1[2] - anisou1[4] * anisou1[4]);
          gradu[2] = gradb * (anisou1[0] * anisou1[1] - anisou1[3] * anisou1[3]);
          gradu[3] = gradb * (2.0 * (anisou1[4] * anisou1[5] - anisou1[3] * anisou1[2]));
          gradu[4] = gradb * (2.0 * (anisou1[3] * anisou1[5] - anisou1[4] * anisou1[1]));
          gradu[5] = gradb * (2.0 * (anisou1[3] * anisou1[4] - anisou1[5] * anisou1[0]));
          a1.addToAnisouGradient(gradu);
          // Atom 2
          double gradBiso = u2b(-gradb * u2 * u2);
          a2.addToTempFactorGradient(gradBiso);
        }
      }
    }
    return e;
  }
}
