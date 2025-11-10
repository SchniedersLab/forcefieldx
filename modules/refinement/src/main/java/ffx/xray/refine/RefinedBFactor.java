package ffx.xray.refine;

import ffx.potential.bonded.Atom;

import java.util.logging.Logger;

import static ffx.numerics.math.MatrixMath.determinant3;
import static ffx.numerics.math.ScalarMath.b2u;
import static ffx.numerics.math.ScalarMath.u2b;
import static java.lang.String.format;

/**
 * The RefinedBfactor class defines an object for handling the refinement
 * of B-factors (either isotropic or anisotropic) for a specific atom.
 * It allows the setting and retrieval of atomic B-factors and anisotropic
 * U-matrix values. Additionally, it supports grouping atoms to maintain
 * constrained B-factor refinements across the group.
 */
public class RefinedBFactor extends RefinedParameter {

  private static final Logger logger = Logger.getLogger(RefinedBFactor.class.getName());

  private static final double BFACTOR_SCALE = 12.0;
  private static final double ANISOU_SCALE = 1200.0;

  /**
   * Indicates whether the atom associated with this instance has anisotropic
   * displacement parameters (ANISOU data) available.
   * <p>
   * If true, the atom has ANISOU data, which provides a more detailed
   * representation of atomic displacement using a 3x3 covariance matrix of atomic
   * positional uncertainty. If false, the atom only has an isotropic displacement
   * parameter (B-factor), which assumes uniform displacement in all directions.
   * <p>
   * This field is final and immutable, being determined at the creation of the
   * RefinableBfactor instance based on the presence of anisotropic data for the
   * corresponding atom.
   */
  private final boolean isAnisou;

  /**
   * Constructor for the RefinableBfactor class.
   * Initializes the RefinableBfactor instance by associating it with the given atom
   * and determining if the atom contains anisotropic B-factor data.
   *
   * @param atom the atom associated with the RefinableBfactor instance. This atom
   *             is checked for anisotropic B-factor data during initialization.
   */
  public RefinedBFactor(Atom atom) {
    super(atom);
    double[] anisou = this.atom.getAnisou(null);
    this.isAnisou = anisou != null;
  }

  /**
   * Adds the specified atom to the list of constrained atoms while applying specific constraints.
   * This method sets the atom to active, assigns its temperature factor based on the associated
   * atom, and, if applicable, copies its anisotropic displacement parameters (ANISOU).
   *
   * @param atom the atom to be constrained and updated with relevant parameters such as the
   *             temperature factor and anisotropic displacement data (if available).
   */
  @Override
  public void addConstrainedAtom(Atom atom) {
    // Apply the constraint.
    atom.setActive(true);
    atom.setTempFactor(this.atom.getTempFactor());
    if (isAnisou) {
      atom.setAnisou(this.atom.getAnisou(null));
    }
    constrainedAtoms.add(atom);
  }

  /**
   * Adds the given atom to the list of constrained atoms that scatter in structural models.
   * Sets the atom as active, assigns its temperature factor based on the associated atom,
   * and if applicable, propagates its anisotropic B-factor data (ANISOU).
   *
   * @param atom the atom to be constrained and included in the scattering group.
   */
  @Override
  public void addConstrainedAtomThatScatters(Atom atom) {
    // Apply the constraint.
    atom.setActive(true);
    atom.setTempFactor(this.atom.getTempFactor());
    if (isAnisou) {
      atom.setAnisou(this.atom.getAnisou(null));
    }
    constrainedAtomsThatScatter.add(atom);
  }

  /**
   * Retrieves the number of parameters associated with this instance,
   * based on whether the atom contains anisotropic B-factor data.
   *
   * @return 6 if anisotropic B-factor data is present; otherwise, 1.
   */
  @Override
  public int getNumberOfParameters() {
    if (isAnisou) {
      return 6;
    }
    return 1;
  }

  /**
   * Determines whether the atom associated with this RefinableBfactor instance contains
   * anisotropic B-factor data.
   *
   * @return true if the atom contains anisotropic B-factor data, false otherwise.
   */
  public boolean isAnisou() {
    return isAnisou;
  }

  /**
   * Sets the B-factor value for the associated atom and all constrained atoms.
   * The B-factor, also referred to as the temperature factor, is a value associated
   * with the atomic displacement or mobility in structural models.
   *
   * @param bFactor the B-factor value to be assigned to the associated atom and all constrained atoms
   */
  public void setBFactor(double bFactor) {
    atom.setTempFactor(bFactor);
    double[] anisou = new double[6];
    if (isAnisou) {
      double u = b2u(bFactor);
      anisou[0] = u;
      anisou[1] = u;
      anisou[2] = u;
      atom.setAnisou(anisou);
    }
    for (Atom a : constrainedAtoms) {
      a.setTempFactor(bFactor);
      if (isAnisou && !a.isHydrogen()) {
        a.setAnisou(anisou);
      }
    }
    for (Atom a : constrainedAtomsThatScatter) {
      a.setTempFactor(bFactor);
      if (isAnisou && !a.isHydrogen()) {
        a.setAnisou(anisou);
      }
    }
  }

  /**
   * Sets the anisotropic displacement parameters (ANISOU) for the associated atom.
   * This method adjusts the temperature factor and ANISOU parameters based on the
   * input values, while ensuring that the determinant of the ANISOU matrix is
   * appropriately handled. If the determinant is negative, default values are used.
   * Additionally, sets the temperature factor and ANISOU for constrained atoms.
   *
   * @param anisou an array of six double values representing the ANISOU matrix parameters.
   *               These values define anisotropic atomic displacement.
   */
  public void setAnisou(double[] anisou) {
    if (!isAnisou) {
      return;
    }
    double det = determinant3(anisou);
    if (det > 0.0) {
      atom.setAnisou(anisou);
      det = Math.pow(det, 1.0 / 3.0);
      atom.setTempFactor(u2b(det));
    } else {
      atom.setTempFactor(1.0);
      double u = b2u(1.0);
      anisou[0] = u;
      anisou[1] = u;
      anisou[2] = u;
      anisou[3] = 0.0;
      anisou[4] = 0.0;
      anisou[5] = 0.0;
      atom.setAnisou(anisou);
      logger.info(" Negative ANISOU: " + atom);
    }
    double bfactor = atom.getTempFactor();
    for (Atom a : constrainedAtoms) {
      a.setTempFactor(bfactor);
      if (!a.isHydrogen()) {
        a.setAnisou(anisou);
      }
    }
    for (Atom a : constrainedAtomsThatScatter) {
      a.setTempFactor(bfactor);
      if (!a.isHydrogen()) {
        a.setAnisou(anisou);
      }
    }
  }

  /**
   * Loads the parameters for the atom based on whether it has anisotropic data (ANISOU) or isotropic
   * temperature factor (B-factor). ANISOU data consists of six parameters, while isotropic data is
   * represented by a single parameter.
   *
   * @param parameters an array of doubles to populate with the relevant parameters.
   */
  @Override
  public void getParameters(double[] parameters) {
    if (isAnisou) {
      double[] anisou = new double[6];
      atom.getAnisou(anisou);
      parameters[index] = anisou[0];
      parameters[index + 1] = anisou[1];
      parameters[index + 2] = anisou[2];
      parameters[index + 3] = anisou[3];
      parameters[index + 4] = anisou[4];
      parameters[index + 5] = anisou[5];
    } else {
      parameters[index] = atom.getTempFactor();
    }
  }

  /**
   * Stores parameter values for the associated atom, determining whether to store
   * anisotropic displacement parameters (ANISOU) or isotropic B-factor values.
   * If the atom contains anisotropic data, six ANISOU parameters are stored;
   * otherwise, a single B-factor value is set.
   *
   * @param parameters an array of doubles containing all refined parameters.
   */
  @Override
  public void setParameters(double[] parameters) {
    if (isAnisou) {
      double[] anisou = new double[6];
      anisou[0] = parameters[index];
      anisou[1] = parameters[index + 1];
      anisou[2] = parameters[index + 2];
      anisou[3] = parameters[index + 3];
      anisou[4] = parameters[index + 4];
      anisou[5] = parameters[index + 5];
      setAnisou(anisou);
    } else {
      setBFactor(parameters[index]);
    }
  }

  /**
   * Populates the provided array with velocity values for the associated atom.
   * Depending on whether the atom has anisotropic B-factor data, either six ANISOU velocity
   * components or a single isotropic B-factor velocity will be stored.
   *
   * @param velocity an array of double values to be updated with velocities.
   *                   If the associated atom has anisotropic B-factor data, six
   *                   velocity values (ANISOU components) are stored. Otherwise,
   *                   a single velocity value (isotropic B-factor velocity) is stored.
   */
  @Override
  public void getVelocity(double[] velocity) {
    if (isAnisou) {
      double[] anisou = new double[6];
      atom.getAnisouVelocity(anisou);
      velocity[index] = anisou[0];
      velocity[index + 1] = anisou[1];
      velocity[index + 2] = anisou[2];
      velocity[index + 3] = anisou[3];
      velocity[index + 4] = anisou[4];
      velocity[index + 5] = anisou[5];
    } else {
      velocity[index] = atom.getTempFactorVelocity();
    }
  }

  /**
   * Stores parameter values for the associated atom, determining whether to store
   * anisotropic displacement parameters (ANISOU) or isotropic B-factor values.
   * If the atom contains anisotropic data, six ANISOU parameters are stored;
   * otherwise, a single B-factor value is set.
   *
   * @param velocity an array of doubles containing all refined parameters.
   */
  @Override
  public void setVelocity(double[] velocity) {
    if (isAnisou) {
      double[] anisou = new double[6];
      anisou[0] = velocity[index];
      anisou[1] = velocity[index + 1];
      anisou[2] = velocity[index + 2];
      anisou[3] = velocity[index + 3];
      anisou[4] = velocity[index + 4];
      anisou[5] = velocity[index + 5];
      atom.setAnisouVelocity(anisou);
    } else {
      atom.setTempFactorVelocity(velocity[index]);
    }
  }

  /**
   * Populates the provided array with acceleration values for the associated atom.
   * If the atom contains anisotropic B-factor data (ANISOU), six acceleration components
   * are stored. Otherwise, a single isotropic acceleration value is stored.
   *
   * @param acceleration an array of double values to be updated with acceleration data.
   *                      If the associated atom has anisotropic B-factor data, six values
   *                      (ANISOU components) are stored. Otherwise, a single isotropic
   *                      acceleration value is stored.
   */
  @Override
  public void getAcceleration(double[] acceleration) {
    if (isAnisou) {
      double[] anisou = new double[6];
      atom.getAnisouAcceleration(anisou);
      acceleration[index] = anisou[0];
      acceleration[index + 1] = anisou[1];
      acceleration[index + 2] = anisou[2];
      acceleration[index + 3] = anisou[3];
      acceleration[index + 4] = anisou[4];
      acceleration[index + 5] = anisou[5];
    } else {
      acceleration[index] = atom.getTempFactorAcceleration();
    }
  }

  /**
   * Sets the acceleration for the atom based on the input array.
   * The behavior changes depending on whether the atom is anisotropic or not.
   *
   * @param acceleration an array of doubles representing acceleration values.
   *                      If the atom is anisotropic, six values are extracted
   *                      to set anisotropic acceleration. Otherwise, a single
   *                      value is used to set the isotropic acceleration.
   */
  @Override
  public void setAcceleration(double[] acceleration) {
    if (isAnisou) {
      double[] anisou = new double[6];
      anisou[0] = acceleration[index];
      anisou[1] = acceleration[index + 1];
      anisou[2] = acceleration[index + 2];
      anisou[3] = acceleration[index + 3];
      anisou[4] = acceleration[index + 4];
      anisou[5] = acceleration[index + 5];
      atom.setAnisouAcceleration(anisou);
    } else {
      atom.setTempFactorAcceleration(acceleration[index]);
    }
  }

  /**
   * Populates the provided array with previous acceleration values for the associated atom.
   * If the atom contains anisotropic B-factor data (ANISOU), six acceleration components
   * are stored. Otherwise, a single isotropic acceleration value is stored.
   *
   * @param previousAcceleration an array of double values to be updated with previous acceleration data.
   *                      If the associated atom has anisotropic B-factor data, six values
   *                      (ANISOU components) are stored. Otherwise, a single isotropic
   *                      acceleration value is stored.
   */
  @Override
  public void getPreviousAcceleration(double[] previousAcceleration) {
    if (isAnisou) {
      double[] anisou = new double[6];
      atom.getAnisouPreviousAcceleration(anisou);
      previousAcceleration[index] = anisou[0];
      previousAcceleration[index + 1] = anisou[1];
      previousAcceleration[index + 2] = anisou[2];
      previousAcceleration[index + 3] = anisou[3];
      previousAcceleration[index + 4] = anisou[4];
      previousAcceleration[index + 5] = anisou[5];
    } else {
      previousAcceleration[index] = atom.getTempFactorPreviousAcceleration();
    }
  }

  /**
   * Sets the previous acceleration for the atom based on the input array.
   * The behavior changes depending on whether the atom is anisotropic or not.
   *
   * @param previousAcceleration an array of doubles representing previous acceleration values.
   *                      If the atom is anisotropic, six values are extracted
   *                      to set anisotropic acceleration. Otherwise, a single
   *                      value is used to set the isotropic acceleration.
   */
  @Override
  public void setPreviousAcceleration(double[] previousAcceleration) {
    if (isAnisou) {
      double[] anisou = new double[6];
      anisou[0] = previousAcceleration[index];
      anisou[1] = previousAcceleration[index + 1];
      anisou[2] = previousAcceleration[index + 2];
      anisou[3] = previousAcceleration[index + 3];
      anisou[4] = previousAcceleration[index + 4];
      anisou[5] = previousAcceleration[index + 5];
      atom.setAnisouPreviousAcceleration(anisou);
    } else {
      atom.setTempFactorPreviousAcceleration(previousAcceleration[index]);
    }
  }

  /**
   * Sets the scaling factors for optimization parameters based on whether the atom contains
   * anisotropic B-factor data. If the atom is identified as having anisotropic displacement
   * data, six scaling factors are set. Otherwise, a single scaling factor is applied.
   *
   * @param optimizationScaling an array of doubles representing the scaling factors for optimization
   *                            parameters. This array will be updated with the appropriate scaling
   *                            values based on the isotropic or anisotropic B-factor data.
   */
  @Override
  public void setOptimizationScaling(double[] optimizationScaling) {
    if (isAnisou) {
      optimizationScaling[index] = ANISOU_SCALE;
      optimizationScaling[index + 1] = ANISOU_SCALE;
      optimizationScaling[index + 2] = ANISOU_SCALE;
      optimizationScaling[index + 3] = ANISOU_SCALE;
      optimizationScaling[index + 4] = ANISOU_SCALE;
      optimizationScaling[index + 5] = ANISOU_SCALE;
    } else {
      optimizationScaling[index] = BFACTOR_SCALE;
    }
  }

  /**
   * The zero the B-factor gradient.
   */
  @Override
  public void zeroGradient() {
    if (!isAnisou) {
      atom.setTempFactorGradient(0.0);
      for (Atom a : constrainedAtoms) {
        a.setTempFactorGradient(0.0);
      }
      for (Atom a : constrainedAtomsThatScatter) {
        a.setTempFactorGradient(0.0);
      }
    } else {
      double[] anisouGrad = new double[6];
      atom.setTempFactorGradient(0.0);
      atom.setAnisouGradient(anisouGrad);
      for (Atom a : constrainedAtoms) {
        a.setTempFactorGradient(0.0);
        if (!a.isHydrogen()) {
          a.setAnisouGradient(anisouGrad);
        }
      }
      for (Atom a : constrainedAtomsThatScatter) {
        a.setTempFactorGradient(0.0);
        if (!a.isHydrogen()) {
          a.setAnisouGradient(anisouGrad);
        }
      }
    }
  }

  /**
   * Adds the gradient contributions for the associated atom and all constrained atoms
   * to the provided gradient array. This method accounts for isotropic or anisotropic
   * B-factor (temperature factor) data depending on the current state of the `isAnisou`
   * field. If the atom has anisotropic B-factor data, the gradient is updated using
   * the corresponding ANISOU components. Otherwise, isotropic B-factor gradients are
   * used.
   *
   * @param gradient an array of doubles representing the gradient values to be updated.
   */
  @Override
  public void getGradient(double[] gradient) {
    if (!isAnisou) {
      double bfactor = atom.getTempFactorGradient();
      gradient[index] = bfactor;
      for (Atom a : constrainedAtomsThatScatter) {
        bfactor = a.getTempFactorGradient();
        gradient[index] += bfactor;
      }
    } else {
      double[] anisou = new double[6];
      atom.getAnisouGradient(anisou);
      gradient[index] = anisou[0];
      gradient[index + 1] = anisou[1];
      gradient[index + 2] = anisou[2];
      gradient[index + 3] = anisou[3];
      gradient[index + 4] = anisou[4];
      gradient[index + 5] = anisou[5];
      for (Atom a : constrainedAtomsThatScatter) {
        if (!a.isHydrogen()) {
          a.getAnisouGradient(anisou);
          gradient[index] += anisou[0];
          gradient[index + 1] += anisou[1];
          gradient[index + 2] += anisou[2];
          gradient[index + 3] += anisou[3];
          gradient[index + 4] += anisou[4];
          gradient[index + 5] += anisou[5];
        } else {
          double bfactor = a.getTempFactorGradient();
          double u = b2u(bfactor);
          gradient[index] += u;
          gradient[index + 1] += u;
          gradient[index + 2] += u;
          // No off-diagonal
        }
      }
    }
  }

  /**
   * Populates the provided array with mass values based on the anisotropic or isotropic
   * nature of the associated atom. If the atom contains anisotropic B-factor data,
   * values are assigned in six positions. Otherwise, only one value is assigned.
   *
   * @param mass an array of double values where the mass information will be stored.
   */
  @Override
  public void getMass(double[] mass, double defaultMass) {
    if (isAnisou) {
      mass[index] = defaultMass;
      mass[index + 1] = defaultMass;
      mass[index + 2] = defaultMass;
      mass[index + 3] = defaultMass;
      mass[index + 4] = defaultMass;
      mass[index + 5] = defaultMass;
    } else {
      mass[index] = defaultMass;
    }
  }

  @Override
  public String toString() {
    if (isAnisou) {
      double[] anisou = new double[6];
      atom.getAnisou(anisou);
      StringBuilder sb = new StringBuilder(" Anisotropic B-factor:");
      for (int i = 0; i < 6; i++) {
        sb.append(format("  %10.6f", anisou[i]));
      }
      return sb.toString();
    } else {
      return format(" B-factor:  %10.6f", atom.getTempFactor());
    }
  }
}
