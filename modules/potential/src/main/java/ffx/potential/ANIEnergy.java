// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.potential;

import com.sun.jna.NativeLong;
import edu.uiowa.torchani.TorchANIUtils;
import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parameters.ForceField;

import java.util.logging.Logger;

import static edu.uiowa.torchani.TorchANILibrary.ctorch;
import static ffx.utilities.Constants.HARTREE_TO_KCAL_PER_MOL;

public class ANIEnergy implements Potential, LambdaInterface {

  /** A Logger for the ANIEnergy class. */
  private static final Logger logger = Logger.getLogger(ANIEnergy.class.getName());

  private final int nAtoms;
  private final Atom[] atoms;
  private final NativeLong nAtomsLong;
  private final int[] species;
  private final double[] coordinates;
  private final double[] grad;
  private final String pathToANI;
  private double lambda = 1.0;
  private double energy;
  private final MolecularAssembly molecularAssembly;

  private STATE state = STATE.FAST;

  public ANIEnergy(MolecularAssembly molecularAssembly) {
    TorchANIUtils.init();
    System.out.println(" ANI Dir: " + TorchANIUtils.getLibDirectory());

    atoms = molecularAssembly.getAtomArray();
    nAtoms = atoms.length;
    nAtomsLong = new NativeLong(nAtoms);
    species = new int[nAtoms];
    coordinates = new double[3 * nAtoms];
    grad = new double[3 * nAtoms];
    int index = 0;
    for (Atom atom : atoms) {
      species[index] = atom.getAtomicNumber();
      index++;
    }

    this.molecularAssembly = molecularAssembly;

    ForceField forceField = molecularAssembly.getForceField();
    pathToANI = forceField.getString("ANI_PATH", "ANI2x.pt");
  }

  /**
   * Compute the ANI energy and gradint.
   *
   * @param gradient If true, compute the gradient.
   * @param print If true, turn on extra printing.
   * @return The ANI energy.
   */
  public double energy(boolean gradient, boolean print) {
    getCoordinates(coordinates);
    if (gradient) {
      energyAndGradient(coordinates, grad);
      for (int i = 0; i < nAtoms; i++) {
        Atom ai = atoms[i];
        int index = 3 * i;
        ai.addToXYZGradient(grad[index], grad[index + 1], grad[index + 2]);
      }
    } else {
      energy(coordinates);
    }
    return energy;
  }

  @Override
  public double energy(double[] x) {
    energy = ctorch(pathToANI, nAtomsLong, species, x, grad) * HARTREE_TO_KCAL_PER_MOL;
    return energy;
  }

  @Override
  public double energyAndGradient(double[] x, double[] g) {
    energy = ctorch(pathToANI, nAtomsLong, species, x, g) * HARTREE_TO_KCAL_PER_MOL;
    for (int i = 0; i < grad.length; i++) {
      grad[i] = g[i] * HARTREE_TO_KCAL_PER_MOL;
    }
    return energy;
  }

  @Override
  public double[] getAcceleration(double[] acceleration) {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    return forceFieldEnergy.getAcceleration(acceleration);
  }

  @Override
  public double[] getCoordinates(double[] parameters) {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    return forceFieldEnergy.getCoordinates(parameters);
  }

  @Override
  public void setCoordinates(double[] parameters) {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    forceFieldEnergy.setCoordinates(parameters);
  }

  @Override
  public STATE getEnergyTermState() {
    return state;
  }

  @Override
  public void setEnergyTermState(STATE state) {
    this.state = state;
  }

  @Override
  public double[] getMass() {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    return forceFieldEnergy.getMass();
  }

  @Override
  public int getNumberOfVariables() {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    return forceFieldEnergy.getNumberOfVariables();
  }

  @Override
  public double[] getPreviousAcceleration(double[] previousAcceleration) {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    return forceFieldEnergy.getPreviousAcceleration(previousAcceleration);
  }

  @Override
  public double[] getScaling() {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    return forceFieldEnergy.getScaling();
  }

  @Override
  public void setScaling(double[] scaling) {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    forceFieldEnergy.setScaling(scaling);
  }

  @Override
  public double getTotalEnergy() {
    return energy;
  }

  @Override
  public VARIABLE_TYPE[] getVariableTypes() {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    return forceFieldEnergy.getVariableTypes();
  }

  @Override
  public double[] getVelocity(double[] velocity) {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    return forceFieldEnergy.getVelocity(velocity);
  }

  @Override
  public void setAcceleration(double[] acceleration) {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    forceFieldEnergy.setAcceleration(acceleration);
  }

  @Override
  public void setPreviousAcceleration(double[] previousAcceleration) {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    forceFieldEnergy.setPreviousAcceleration(previousAcceleration);
  }

  @Override
  public void setVelocity(double[] velocity) {
    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    forceFieldEnergy.setVelocity(velocity);
  }

  @Override
  public double getLambda() {
    return lambda;
  }

  @Override
  public void setLambda(double lambda) {
    this.lambda = lambda;
  }

  @Override
  public double getd2EdL2() {
    return 0;
  }

  @Override
  public double getdEdL() {
    return energy;
  }

  @Override
  public void getdEdXdL(double[] gradient) {
    System.arraycopy(this.grad, 0, gradient, 0, gradient.length);
  }

}
