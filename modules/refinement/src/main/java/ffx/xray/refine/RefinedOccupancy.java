package ffx.xray.refine;

import ffx.potential.bonded.Atom;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.asin;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

/**
 * The RefinableOccupancy class represents an atom and manages its occupancy value
 * along with any constrained atoms that share the same occupancy constraints.
 * It allows setting and retrieving the occupancy value, ensuring that all
 * constrained atoms are synchronized with the primary atom's occupancy.
 */
public class RefinedOccupancy extends RefinedParameter {

  private static final double OCCUPANCY_SCALE = 120.0;

  /**
   * A collection of atoms that are constrained in relation to a primary {@link Atom}.
   * These constrained atoms are defined as part of the experimental refinement process
   * and are related to the primary atom based on specific refinement requirements.
   * <p>
   * The refined parameter(s) will be applied to the constrained Atoms, but because these
   * atoms do not contribute to experimental scattering, their partial derivatives do not
   * need to be considered.
   * <p>
   * The occupancy of these atoms will be set to 1.0 - occupancy of this parameter.
   */
  protected final List<Atom> complementList;

  /**
   * A collection of {@link Atom} instances that are both constrained with respect to a primary {@link Atom}
   * and directly contribute to experimental scattering.
   * <p>
   * The partial derivative of the experimental refinement target with respect to this
   * parameter(s) is contributed to by each constrained atom that scatters.
   * <p>
   * The occupancy of these atoms will be set to 1.0 - occupancy of this parameter.
   * There contribution to the gradient will be negated.
   */
  protected final List<Atom> complementScatterList;

  /**
   * Constructs a RefinableOccupancy instance for managing the specified atom.
   *
   * @param atom the primary atom whose occupancy value will be managed and potentially constrained
   *             to synchronize with other atoms
   */
  public RefinedOccupancy(Atom atom) {
    super(atom);
    complementList = new ArrayList<>();
    complementScatterList = new ArrayList<>();
  }

  @Override
  public void addConstrainedAtom(Atom atom) {
    // Apply the constraint.
    atom.setActive(true);
    atom.setOccupancy(getOccupancy());
    constrainedAtoms.add(atom);
  }

  @Override
  public void addConstrainedAtomThatScatters(Atom atom) {
    // Apply the constraint.
    atom.setActive(true);
    atom.setOccupancy(getOccupancy());
    constrainedAtomsThatScatter.add(atom);
  }

  /**
   * Adds an atom constrained to this RefinedOccupancy instance. The specified
   * atom will have its active status enabled and its occupancy value set to complement the
   * occupancy of the primary atom managed by this instance.
   *
   * @param atom the atom to be added as a complementary, constrained atom. Its occupancy will
   *             be adjusted to ensure it complements the occupancy of the primary atom.
   */
  public void addConstrainedAtomComplement(Atom atom) {
    // Apply the constraint.
    atom.setActive(true);
    atom.setOccupancy(1.0 - getOccupancy());
    complementList.add(atom);
  }

  /**
   * Adds an atom constrained to this RefinedOccupancy instance as a complementary, scattering atom.
   * The specified atom will have its active status enabled and its occupancy value set to complement
   * the occupancy of the primary atom managed by this instance. Additionally, the atom will be added
   * to the complementScatterList.
   *
   * @param atom the atom to be added as a complementary, constrained, scattering atom. Its occupancy
   *             will be adjusted to ensure it complements the occupancy of the primary atom.
   */
  public void addConstrainedAtomThatScattersComplement(Atom atom) {
    // Apply the constraint.
    atom.setActive(true);
    atom.setOccupancy(1.0 - getOccupancy());
    complementScatterList.add(atom);
  }

  @Override
  public int getNumberOfParameters() {
    return 1;
  }

  /**
   * Retrieves the occupancy value of the primary atom managed by this instance.
   * The occupancy represents the fractional occupancy of the atom in the molecular structure.
   *
   * @return the occupancy value of the primary atom as a double
   */
  public double getOccupancy() {
    return atom.getOccupancy();
  }

  /**
   * Sets the occupancy value for the primary atom and any constrained atoms, ensuring
   * that all atoms share the same occupancy value. The occupancy represents the
   * fractional occupancy of the atom in the molecular structure.
   *
   * @param occupancy the new occupancy value to be assigned to the primary atom and
   *                  all constrained atoms
   */
  public void setOccupancy(double occupancy) {
    atom.setOccupancy(occupancy);
    for (Atom a : constrainedAtoms) {
      a.setOccupancy(occupancy);
    }
    for (Atom a : constrainedAtomsThatScatter) {
      a.setOccupancy(occupancy);
    }
    for (Atom a : complementList) {
      a.setOccupancy(1.0 - occupancy);
    }
    for (Atom a : complementScatterList) {
      a.setOccupancy(1.0 - occupancy);
    }
  }

  /**
   * Updates the specified parameters array with the occupancy value of the primary atom
   * managed by this instance at the specified index.
   *
   * @param parameters an array of doubles representing the parameters to be updated.
   */
  @Override
  public void getParameters(double[] parameters) {
    double occupancy = atom.getOccupancy();
    // Convert occupancy to theta:
    // occupancy = sin^2(theta)
    // theta = asin(sqrt(occupancy);
    parameters[index] = asin(sqrt(occupancy));
  }

  /**
   * Sets the occupancy value for the primary atom and any
   * constrained atoms using a specified index within the parameters array.
   *
   * @param parameters an array of doubles containing all refined parameters.
   */
  @Override
  public void setParameters(double[] parameters) {
    // occupancy = sin^2(theta)
    double sinTheta = sin(parameters[index]);
    double occupancy = sinTheta * sinTheta;
    setOccupancy(occupancy);
  }

  /**
   * Updates the specified parameters array with the occupancy velocity of the primary atom
   * managed by this instance at the specified index. The velocity represents the rate of change
   * of the atom's occupancy with units of degrees per picosecond.
   *
   * @param velocity an array of doubles to be updated with the atom's occupancy velocity.
   */
  @Override
  public void getVelocity(double[] velocity) {
    // Occupancy is stored from 0 to 1: occupancy = sin^2(theta).
    // Theta is sent to the integrator / optimizer: theta = asin(sqrt(occupancy)).
    // The occupancy "particle" has velocity with units of degrees / picosecond.
    double v = atom.getOccupancyVelocity();
    velocity[index] = v;
  }

  /**
   * Updates the occupancy velocity for the primary atom managed by this instance
   * using a specified value from the parameters array at the specified index.
   *
   * @param velocity an array of doubles representing refined parameters,
   *                 where the relevant velocity value is extracted and applied.
   */
  @Override
  public void setVelocity(double[] velocity) {
    atom.setOccupancyVelocity(velocity[index]);
  }

  /**
   * Updates the specified acceleration array with the occupancy acceleration
   * of the primary atom managed by this instance at the relevant index.
   * The acceleration represents the second derivative of the atom's occupancy
   * with respect to time.
   *
   * @param acceleration an array of doubles to be updated with the atom's occupancy acceleration.
   */
  @Override
  public void getAcceleration(double[] acceleration) {
    double v = atom.getOccupancyAcceleration();
    acceleration[index] = v;
  }

  /**
   * Sets the occupancy acceleration for the primary atom based on the specified array.
   * This method updates the occupancy acceleration of the atom managed by this instance
   * using the value at the specified index.
   *
   * @param acceleration an array of doubles representing the acceleration values
   *                     to be used for updating the atom's occupancy acceleration.
   */
  @Override
  public void setAcceleration(double[] acceleration) {
    atom.setOccupancyAcceleration(acceleration[index]);
  }

  /**
   * Updates the specified previous acceleration array with the occupancy acceleration
   * of the primary atom managed by this instance at the relevant index.
   * The previous acceleration represents the second derivative of the atom's occupancy
   * with respect to time.
   *
   * @param previousAcceleration an array of doubles to be updated with the atom's occupancy acceleration.
   */
  @Override
  public void getPreviousAcceleration(double[] previousAcceleration) {
    double v = atom.getOccupancyAcceleration();
    previousAcceleration[index] = v;
  }

  /**
   * Sets the occupancy previous acceleration for the primary atom based on the specified array.
   * This method updates the occupancy previous acceleration of the atom managed by this instance
   * using the value at the specified index.
   *
   * @param previousAcceleration an array of doubles representing the previous acceleration values
   *                             to be used for updating the atom's occupancy acceleration.
   */
  @Override
  public void setPreviousAcceleration(double[] previousAcceleration) {
    atom.setOccupancyAcceleration(previousAcceleration[index]);
  }

  @Override
  public void setOptimizationScaling(double[] optimizationScaling) {
    optimizationScaling[index] = OCCUPANCY_SCALE;
  }

  /**
   * Resets the occupancy gradient associated with the primary atom and any constrained
   * atoms to zero. This method ensures that the occupancy gradient for the primary
   * atom, constrained atoms, and constrained atoms that scatter are all reset to a value of 0.0.
   * <p>
   * The occupancy gradient represents the partial derivative of the energy with
   * respect to the atom's occupancy, and this reset is part of the gradient update process.
   */
  @Override
  public void zeroGradient() {
    atom.setOccupancyGradient(0.0);
    for (Atom a : constrainedAtoms) {
      a.setOccupancyGradient(0.0);
    }
    for (Atom a : constrainedAtomsThatScatter) {
      a.setOccupancyGradient(0.0);
    }
    for (Atom a : complementList) {
      a.setOccupancyGradient(0.0);
    }
    for (Atom a : complementScatterList) {
      a.setOccupancyGradient(0.0);
    }
  }

  /**
   * Updates the specified gradient array with the occupancy gradients of the primary atom and any
   * constrained atoms that scatter.
   *
   * @param gradient an array of doubles representing the gradient to be updated.
   */
  @Override
  public void getGradient(double[] gradient) {
    gradient[index] = atom.getOccupancyGradient();
    for (Atom a : constrainedAtomsThatScatter) {
      gradient[index] += a.getOccupancyGradient();
    }
    for (Atom a : complementScatterList) {
      gradient[index] -= a.getOccupancyGradient();
    }
    // Add the chain rule term:
    // dE / dTheta = dE / dOcc * dOcc / dTheta
    // dE / dTheta = dE / dOcc * sin(2 * theta)
    double theta = asin(sqrt(getOccupancy()));
    double sin2Theta = sin(2.0 * theta);
    gradient[index] *= sin2Theta;
  }

  /**
   * Updates the specified mass array with the default mass value at the relevant index.
   *
   * @param mass        an array of doubles representing the mass values to be updated.
   * @param defaultMass a double representing the default mass value to assign.
   */
  @Override
  public void getMass(double[] mass, double defaultMass) {
    mass[index] = defaultMass;
  }

  @Override
  public String toString() {
    return format(" Occupancy: %10.6f", atom.getOccupancy());
  }
}
