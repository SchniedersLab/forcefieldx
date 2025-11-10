package ffx.xray.refine;

import ffx.potential.bonded.Atom;

import javax.annotation.Nullable;

import static java.lang.String.format;

/**
 * Represents a set of coordinates that can be refined during a refinement process.
 * This class encapsulates an {@link Atom} object and any additional atoms that are constrained
 * to follow its coordinates.
 */
public class RefinedCoordinates extends RefinedParameter {

  /**
   * Optimization scale factor for atomic coordinates.
   */
  private static final double COORDINATE_SCALE = 12.0;

  /**
   * Constructs a new RefinableCoordinates instance for the specified atom.
   * This instance represents a set of coordinates that can be refined,
   * with the given atom serving as the primary reference point during the refinement process.
   * Changes to the primary atom's coordinates will also affect any constrained atoms
   * associated with this instance.
   *
   * @param atom the {@link Atom} whose coordinates are being refined; this atom serves
   *             as the primary reference point for the refinement process.
   */
  public RefinedCoordinates(Atom atom) {
    super(atom);
  }

  @Override
  public void addConstrainedAtom(Atom atom) {
    // Apply the constraint.
    atom.setActive(true);
    atom.setXYZ(getCoordinates(null));
    constrainedAtoms.add(atom);
  }

  @Override
  public void addConstrainedAtomThatScatters(Atom atom) {
    // Apply the constraint.
    atom.setActive(true);
    atom.setXYZ(getCoordinates(null));
    constrainedAtomsThatScatter.add(atom);
  }

  /**
   * Retrieves the XYZ coordinates of the primary atom, updating the provided array or
   * returning a new array if the input array is null.
   *
   * @param xyz an array of doubles where the XYZ coordinates should be stored. If the input
   *            array is null, a new double array will be created and returned.
   * @return an array of doubles containing the XYZ coordinates of the primary atom. If the input
   * array is not null, it will be updated with the coordinates and returned.
   */
  public double[] getCoordinates(@Nullable double[] xyz) {
    return atom.getXYZ(xyz);
  }

  /**
   * Sets the coordinates of the primary atom and updates the coordinates of
   * all constrained atoms to match the new values.
   *
   * @param xyz an array of doubles representing the XYZ coordinates to be set for the primary atom and its constrained atoms
   */
  public void setCoordinates(double[] xyz) {
    atom.setXYZ(xyz);
    for (Atom a : constrainedAtoms) {
      a.setXYZ(xyz);
    }
    for (Atom a : constrainedAtomsThatScatter) {
      a.setXYZ(xyz);
    }
  }

  @Override
  public int getNumberOfParameters() {
    return 3;
  }

  /**
   * Updates a subset of the parameter array with the XYZ coordinates of the primary atom.
   *
   * @param parameters an array of doubles where the XYZ coordinates of the primary atom
   *                   will be stored. The x, y, and z values are assigned sequentially
   *                   starting at the current index value.
   */
  @Override
  public void getParameters(double[] parameters) {
    parameters[index] = atom.getX();
    parameters[index + 1] = atom.getY();
    parameters[index + 2] = atom.getZ();
  }

  /**
   * Updates the XYZ coordinates of the primary atom and all associated constrained atoms
   * using the provided parameters array.
   *
   * @param parameters an array of doubles containing all refined parameters.
   */
  @Override
  public void setParameters(double[] parameters) {
    double[] xyz = new double[3];
    xyz[0] = parameters[index];
    xyz[1] = parameters[index + 1];
    xyz[2] = parameters[index + 2];
    setCoordinates(xyz);
  }

  /**
   * Updates a subset of the parameter array with the XYZ velocity components of the primary atom.
   *
   * @param velocity an array of doubles where the XYZ velocity components of the primary atom
   *                 will be stored. The x, y, and z components are assigned sequentially
   *                 starting at the current index value.
   */
  @Override
  public void getVelocity(double[] velocity) {
    double[] v = new double[3];
    atom.getVelocity(v);
    velocity[index] = v[0];
    velocity[index + 1] = v[1];
    velocity[index + 2] = v[2];
  }

  /**
   * Sets the velocity components of the primary atom using the provided parameters array.
   * The velocity is updated with values sequentially extracted from the array at the current index.
   *
   * @param velocity an array of doubles containing velocity data. The x, y, and z components
   *                 of the velocity are taken sequentially starting at the current index.
   */
  @Override
  public void setVelocity(double[] velocity) {
    double[] v = new double[3];
    v[0] = velocity[index];
    v[1] = velocity[index + 1];
    v[2] = velocity[index + 2];
    atom.setVelocity(v);
  }

  /**
   * Updates a subset of the parameter array with the XYZ acceleration components of the primary atom.
   *
   * @param acceleration an array of doubles where the XYZ acceleration components of the primary atom
   *                     will be stored. The x, y, and z components are assigned sequentially
   *                     starting at the current index value.
   */
  @Override
  public void getAcceleration(double[] acceleration) {
    double[] a = new double[3];
    atom.getAcceleration(a);
    acceleration[index] = a[0];
    acceleration[index + 1] = a[1];
    acceleration[index + 2] = a[2];
  }

  /**
   * Sets the acceleration components for the primary atom.
   * The method extracts the x, y, and z components of acceleration
   * from the provided array starting at a specified index and
   * assigns them to the primary atom.
   *
   * @param acceleration an array of doubles containing the XYZ acceleration
   *                     components. The x, y, and z components are retrieved
   *                     sequentially from the array starting at the current index.
   */
  @Override
  public void setAcceleration(double[] acceleration) {
    double[] a = new double[3];
    a[0] = acceleration[index];
    a[1] = acceleration[index + 1];
    a[2] = acceleration[index + 2];
    atom.setAcceleration(a);
  }

  /**
   * Updates a subset of the parameter array with the XYZ previous acceleration components of the primary atom.
   *
   * @param previousAcceleration an array of doubles where the XYZ previous acceleration components
   *                             of the primary atom will be stored. The x, y, and z components
   *                             are assigned sequentially starting at the current index value.
   */
  @Override
  public void getPreviousAcceleration(double[] previousAcceleration) {
    double[] a = new double[3];
    atom.getPreviousAcceleration(a);
    previousAcceleration[index] = a[0];
    previousAcceleration[index + 1] = a[1];
    previousAcceleration[index + 2] = a[2];
  }

  /**
   * Sets the previous acceleration components for the primary atom.
   * The method extracts the x, y, and z components of previous acceleration
   * from the provided array starting at a specified index and
   * assigns them to the primary atom.
   *
   * @param previousAcceleration an array of doubles containing the XYZ previous acceleration
   *                             components. The x, y, and z components are retrieved
   *                             sequentially from the array starting at the current index.
   */
  @Override
  public void setPreviousAcceleration(double[] previousAcceleration) {
    double[] a = new double[3];
    a[0] = previousAcceleration[index];
    a[1] = previousAcceleration[index + 1];
    a[2] = previousAcceleration[index + 2];
    atom.setPreviousAcceleration(a);
  }

  /**
   * Sets the scaling factors for optimization parameters at specific indices.
   * This method updates a portion of the provided array with predefined scaling values.
   *
   * @param optimizationScaling an array of doubles representing scaling factors for optimization.
   *                            The method modifies specific indices within this array.
   */
  @Override
  public void setOptimizationScaling(double[] optimizationScaling) {
    optimizationScaling[index] = COORDINATE_SCALE;
    optimizationScaling[index + 1] = COORDINATE_SCALE;
    optimizationScaling[index + 2] = COORDINATE_SCALE;
  }

  /**
   * Initialize the coordinate gradient to zero.
   */
  @Override
  public void zeroGradient() {
    atom.setXYZGradient(0.0, 0.0, 0.0);
    atom.setLambdaXYZGradient(0.0, 0.0, 0.0);
    for (Atom a : constrainedAtoms) {
      a.setXYZGradient(0.0, 0.0, 0.0);
      a.setLambdaXYZGradient(0.0, 0.0, 0.0);
    }
    for (Atom a : constrainedAtomsThatScatter) {
      a.setXYZGradient(0.0, 0.0, 0.0);
      a.setLambdaXYZGradient(0.0, 0.0, 0.0);
    }
  }

  /**
   * Updates the gradient array with the contributions from the primary atom
   * and any constrained atoms that scatter. The XYZ gradient values are
   * retrieved for each relevant atom and added to the corresponding indices
   * of the provided gradient array.
   *
   * @param gradient an array of doubles where the gradient contributions
   *                 will be accumulated.
   */
  @Override
  public void getGradient(double[] gradient) {
    double[] xyz = new double[3];
    atom.getXYZGradient(xyz);
    gradient[index] = xyz[0];
    gradient[index + 1] = xyz[1];
    gradient[index + 2] = xyz[2];
    for (Atom a : constrainedAtomsThatScatter) {
      a.getXYZGradient(xyz);
      gradient[index] += xyz[0];
      gradient[index + 1] += xyz[1];
      gradient[index + 2] += xyz[2];
    }
  }

  /**
   * Updates the provided array with the atomic mass of the primary atom.
   *
   * @param mass        an array of doubles where the atomic mass will be stored.
   * @param defaultMass the default mass is ignored for coordinates.
   */
  @Override
  public void getMass(double[] mass, double defaultMass) {
    double m = atom.getMass();
    mass[index] = m;
    mass[index + 1] = m;
    mass[index + 2] = m;
  }

  @Override
  public String toString() {
    double[] xyz = new double[3];
    atom.getXYZ(xyz);
    StringBuilder sb = new StringBuilder(" Atomic Coordinates: ");
    for (int i = 0; i < 3; i++) {
      sb.append(format("  %10.6f", xyz[i]));
    }
    return sb.toString();
  }
}
