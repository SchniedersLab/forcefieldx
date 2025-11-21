package ffx.xray.refine;

import ffx.potential.bonded.Atom;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents a parameter that can be refined in an experimental refinement process.
 * Encapsulates a primary {@link Atom} and a collection of constrained atoms that are
 * related to this primary atom based on the specific refinement requirement.
 * The concrete behavior and details of the refinement process are defined by subclasses.
 */
public abstract class RefinedParameter {

  /**
   * The primary {@link Atom} that serves as the central reference in a refinement process.
   * This atom is the focal point around which constrained atoms are defined or managed.
   */
  protected final Atom atom;

  /**
   * The index of this parameter in the overall parameter array. If this parameter is
   * described by more than one variable, the index is for the first variable.
   */
  protected int index;

  /**
   * A collection of atoms that are constrained in relation to a primary {@link Atom}.
   * These constrained atoms are defined as part of the experimental refinement process
   * and are related to the primary atom based on specific refinement requirements.
   * <p>
   * The refined parameter(s) will be applied to the constrained Atoms, but because these
   * atoms do not contribute to experimental scattering, their partial derivatives do not
   * need to be considered.
   */
  protected final List<Atom> constrainedAtoms;

  /**
   * A collection of {@link Atom} instances that are both constrained with respect to a primary {@link Atom}
   * and directly contribute to experimental scattering.
   * <p>
   * The partial derivative of the experimental refinement target with respect to this
   * parameter(s) is contributed to by each constrained atom that scatters.
   */
  protected final List<Atom> constrainedAtomsThatScatter;

  /**
   * Constructs a new RefinedParameters object for the specified primary {@link Atom}.
   * The primary atom is set to active, and initialized lists are created for constrained
   * atoms and constrained atoms that contribute to scattering in the refinement process.
   *
   * @param atom the primary {@link Atom} that serves as the central reference in the refinement process.
   */
  public RefinedParameter(Atom atom) {
    this.atom = atom;
    this.atom.setActive(true);
    this.constrainedAtoms = new ArrayList<>();
    this.constrainedAtomsThatScatter = new ArrayList<>();
  }

  /**
   * Retrieves the primary atom associated with this RefinedParameter instance.
   * The returned atom serves as the central reference in the refinement process.
   *
   * @return the primary {@link Atom} associated with this RefinedParameter.
   */
  public Atom getAtom() {
    return atom;
  }

  /**
   * Sets the index for this RefinedParameter instance.
   *
   * @param index the integer value to set as the index
   */
  public void setIndex(int index) {
    this.index = index;
  }

  /**
   * Retrieves the index value associated with this RefinedParameter instance.
   *
   * @return the integer index value of this RefinedParameter
   */
  public int getIndex() {
    return index;
  }

  /**
   * Adds an atom to the list of constrained atoms for the given primary atom.
   * This method defines the relationship between the primary atom and other atoms
   * that are constrained or related as part of the refinement process.
   *
   * @param atom the atom to be added as a constrained atom associated with the primary atom.
   */
  public abstract void addConstrainedAtom(Atom atom);

  /**
   * Adds an atom to the list of constrained atoms that contribute to experimental scattering.
   *
   * @param atom the atom to be added as a constrained atom that has scattering contributions.
   */
  public abstract void addConstrainedAtomThatScatters(Atom atom);

  /**
   * Retrieves the number of parameters associated with this refineable.
   * Concrete implementations define how the parameters are determined and counted.
   *
   * @return the number of parameters related to this refinement process.
   */
  public abstract int getNumberOfParameters();

  /**
   * Get parameters from this RefinedParameter instance.
   * Concrete implementations define how the given parameter values are processed
   * and associated with the refinement process.
   *
   * @param parameters an array of double values representing all refined parameters.
   */
  public abstract void getParameters(double[] parameters);

  /**
   * Set parameter values into this RefinedParameter instance.
   *
   * @param parameters an array of double values representing all refined parameters.
   */
  public abstract void setParameters(double[] parameters);

  /**
   * Retrieves the velocities associated with the refined parameters.
   *
   * @param velocity an array of double values where the velocities of
   *                   the refined parameters will be stored.
   */
  public abstract void getVelocity(double[] velocity);

  /**
   * Sets the velocities for the refined parameters.
   *
   * @param velocity an array of double values representing the velocities for each refined parameter
   */
  public abstract void setVelocity(double[] velocity);

  /**
   * Retrieves the acceleration values associated with the refined parameters.
   * Concrete implementations define how the acceleration data is retrieved and processed.
   *
   * @param acceleration an array of double values where the accelerations of
   *                     the refined parameters will be stored
   */
  public abstract void getAcceleration(double[] acceleration);

  /**
   * Sets the acceleration values for the refined parameters.
   * This method updates the acceleration for each refined parameter based on the provided array.
   *
   * @param acceleration an array of double values representing the acceleration for each refined parameter
   */
  public abstract void setAcceleration(double[] acceleration);

  /**
   * Retrieves the previous acceleration values associated with the refined parameters.
   * Concrete implementations define how the previous acceleration data is retrieved and processed.
   *
   * @param acceleration an array of double values where the previous accelerations of
   *                     the refined parameters will be stored
   */
  public abstract void getPreviousAcceleration(double[] acceleration);

  /**
   * Sets the previous acceleration values for the refined parameters.
   * This method updates the previous acceleration for each refined parameter based on the provided array.
   *
   * @param acceleration an array of double values representing the previous acceleration for each refined parameter
   */
  public abstract void setPreviousAcceleration(double[] acceleration);

  /**
   * Sets the optimization scaling factors for this RefinedParameter instance.
   * The scaling factors are typically used to adjust the magnitude of parameter
   * updates during optimization processes, ensuring appropriate convergence behavior.
   *
   * @param optimizationScaling an array of double values representing the scaling
   *                            factors for each optimization parameter.
   */
  public abstract void setOptimizationScaling(double[] optimizationScaling);

  /**
   * Zero out the gradient for this RefinedParameter.
   */
  public abstract void zeroGradient();

  /**
   * Load the gradient for this RefinedParameter.
   */
  public abstract void getGradient(double[] gradient);

  /**
   * Mass for extended Lagrangian
   *
   * @param mass Store the mass to use for this parameter.
   * @param defaultMass The default mass to use for this parameter.
   */
  public abstract void getMass(double[] mass, double defaultMass);
}
