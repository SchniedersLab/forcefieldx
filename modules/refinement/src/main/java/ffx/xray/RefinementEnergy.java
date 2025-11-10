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
package ffx.xray;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.thermostats.Thermostat;
import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.realspace.RealSpaceData;
import ffx.realspace.RealSpaceEnergy;
import ffx.xray.refine.RefinementMode;
import ffx.xray.refine.RefinementModel;

import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.kB;
import static java.lang.String.format;
import static java.util.Arrays.fill;

/**
 * Combine the X-ray target and chemical potential energy using the {@link
 * ffx.crystal.CrystalPotential} interface
 *
 * @author Timothy D. Fenn
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class RefinementEnergy implements LambdaInterface, CrystalPotential, AlgorithmListener {

  private static final Logger logger = Logger.getLogger(RefinementEnergy.class.getName());

  /**
   * Container to store experimental data.
   */
  private final DataContainer data;
  /**
   * Refinement model to use.
   */
  private final RefinementModel refinementModel;
  /**
   * Specifies the refinement mode used in the energy refinement process.
   * Controls the approach or strategy applied during optimization or
   * analysis of refinement parameters.
   */
  private final RefinementMode refinementMode;
  /**
   * Total potential energy.
   */
  private double totalEnergy;
  /**
   * An array of atoms being refined.
   */
  private final Atom[] scatteringAtoms;
  /**
   * The number of atoms being refined.
   */
  private final int nAtoms;
  /**
   * MolecularAssembly instances being refined.
   */
  private final MolecularAssembly[] molecularAssemblies;
  /**
   * Atomic coordinates for computing the chemical energy.
   */
  private final double[][] xChemical;
  /**
   * Array for storing chemical gradient.
   */
  private final double[][] gChemical;
  /**
   * The Potential based on experimental data.
   */
  private CrystalPotential dataEnergy;
  /**
   * Array for storing the experimental gradient.
   */
  private double[] gExperiment;
  /**
   * Optimization scale factors.
   */
  private double[] optimizationScaling;
  /**
   * The total number of parameters being refined.
   */
  private final int n;
  /**
   * The number of XYZ coordinates being refined.
   */
  private final int nXYZ;
  /**
   * The number of b-factor parameters being refined.
   */
  private final int nBFactor;
  /**
   * The number of occupancy parameters being refined.
   */
  private final int nOccupancy;
  /**
   * Compute fast varying forces, slowly varying forces, or both.
   */
  private STATE state = STATE.BOTH;

  /**
   * A thermostat instance.
   */
  protected Thermostat thermostat;
  /**
   * The kT scale factor.
   */
  private double kTScale;

  /**
   * Constructs a RefinementEnergy instance with the input data and optimization scaling factors.
   *
   * @param dataContainer the data container that provides refinement-related data such as
   *                      refinement model, scattering atoms, and molecular assemblies
   */
  public RefinementEnergy(DataContainer dataContainer) {
    this.data = dataContainer;
    refinementModel = dataContainer.getRefinementModel();
    refinementMode = refinementModel.getRefinementMode();
    molecularAssemblies = refinementModel.getMolecularAssemblies();
    scatteringAtoms = refinementModel.getScatteringAtoms();
    nAtoms = scatteringAtoms.length;

    thermostat = null;
    kTScale = 1.0;

    // Set the number of refinement parameters.
    nXYZ = refinementModel.getNumCoordParameters();
    nBFactor = refinementModel.getNumBFactorParameters();
    nOccupancy = refinementModel.getNumOccupancyParameters();
    n = nXYZ + nBFactor + nOccupancy;

    // Initialize force field and Xray energies
    for (MolecularAssembly molecularAssembly : molecularAssemblies) {
      ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
      if (forceFieldEnergy == null) {
        forceFieldEnergy = ForceFieldEnergy.energyFactory(molecularAssembly);
        molecularAssembly.setPotential(forceFieldEnergy);
      }
    }

    if (dataContainer instanceof DiffractionData diffractionData) {
      if (!diffractionData.getScaled()[0]) {
        diffractionData.printStats();
      }
      dataEnergy = new XRayEnergy(diffractionData);
      // We will handle parameter (un)scaling within this class.
      dataEnergy.setScaling(null);
    } else if (dataContainer instanceof RealSpaceData realSpaceData) {
      dataEnergy = new RealSpaceEnergy(realSpaceData);
      // We will handle parameter (un)scaling within this class.
      dataEnergy.setScaling(null);
    }

    int assemblySize = molecularAssemblies.length;
    xChemical = new double[assemblySize][];
    gChemical = new double[assemblySize][];
    for (int i = 0; i < assemblySize; i++) {
      int len = molecularAssemblies[i].getActiveAtomArray().length * 3;
      xChemical[i] = new double[len];
      gChemical[i] = new double[len];
    }
    gExperiment = new double[n];
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean algorithmUpdate(MolecularAssembly active) {
    if (thermostat != null) {
      kTScale = KCAL_TO_GRAM_ANG2_PER_PS2 / (thermostat.getTargetTemperature() * kB);
    }
    logger.info(" kTscale: " + kTScale);
    logger.info(data.printEnergyUpdate());
    return true;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean destroy() {
    return dataEnergy.destroy();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energy(double[] x) {
    return energy(x, false);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energy(double[] x, boolean print) {
    double weight = data.getWeight();
    double e = 0.0;

    if (thermostat != null) {
      kTScale = KCAL_TO_GRAM_ANG2_PER_PS2 / (thermostat.getTargetTemperature() * kB);
    }

    unscaleCoordinates(x);
    refinementModel.setParameters(x);
    RefinementMode refinementMode = refinementModel.getRefinementMode();

    int numAssemblies = molecularAssemblies.length;

    if (refinementMode.includesCoordinates()) {
      // Compute the chemical energy.
      for (int i = 0; i < numAssemblies; i++) {
        ForceFieldEnergy forceFieldEnergy = molecularAssemblies[i].getPotentialEnergy();
        forceFieldEnergy.getCoordinates(xChemical[i]);
        double curE = forceFieldEnergy.energy(xChemical[i], print);
        e += curE;
      }
      e = e * kTScale / numAssemblies;
      // Compute the experimental target energy.
      e += weight * dataEnergy.energy(x, print);
    } else {
      // Only compute the experimental target energy.
      e = dataEnergy.energy(x, print);
    }

    scaleCoordinates(x);

    totalEnergy = e;
    return e;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Implementation of the {@link CrystalPotential} interface for the RefinementEnergy.
   */
  @Override
  public double energyAndGradient(double[] x, double[] g) {
    return energyAndGradient(x, g, false);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Implementation of the {@link CrystalPotential} interface for the RefinementEnergy.
   */
  @Override
  public double energyAndGradient(double[] x, double[] g, boolean print) {
    double weight = data.getWeight();
    double e = 0.0;
    fill(g, 0.0);
    fill(gExperiment, 0.0);

    if (thermostat != null) {
      kTScale = KCAL_TO_GRAM_ANG2_PER_PS2 / (thermostat.getTargetTemperature() * kB);
    }

    unscaleCoordinates(x);
    refinementModel.setParameters(x);

    if (refinementMode.includesCoordinates()) {
      int numAssemblies = molecularAssemblies.length;
      // Compute the chemical energy and gradient.
      for (int i = 0; i < numAssemblies; i++) {
        ForceFieldEnergy forceFieldEnergy = molecularAssemblies[i].getPotentialEnergy();
        forceFieldEnergy.getCoordinates(xChemical[i]);
        double curE = forceFieldEnergy.energyAndGradient(xChemical[i], gChemical[i], print);
        e += curE;
        // Aggregate the contribution of this conformer into the overall gradient.
        refinementModel.addAssemblyGradient(i, g);
      }

      e = kTScale * e / numAssemblies;
      // normalize gradient for multiple-counted atoms
      if (numAssemblies > 1) {
        for (int i = 0; i < nXYZ; i++) {
          g[i] /= numAssemblies;
        }
      }
      for (int i = 0; i < nXYZ; i++) {
        g[i] *= kTScale;
      }

      double xE = dataEnergy.energyAndGradient(x, gExperiment);
      e += weight * xE;

      // Add the chemical coordinate gradient to experimental coordinate gradient.
      for (int i = 0; i < nXYZ; i++) {
        g[i] += weight * gExperiment[i];
      }

//      for (int i = 0; i < 10; i++) {
//        Atom atom = scatteringAtoms[i];
//        double width = atom.getFormFactorWidth();
//        double bfactor = atom.getTempFactor();
//        double[] grad = new double[3];
//        atom.getXYZGradient(grad);
//        logger.info(format(" Atom Grad %d w=%16.8f b=%16.8f: %16.8f %16.8f %16.8f",
//            i, width, bfactor, grad[0], grad[1], grad[2]));
//        grad[0] = gExperiment[i * 3];
//        grad[1] = gExperiment[i * 3 + 1];
//        grad[2] = gExperiment[i * 3 + 2];
//        logger.info(format(" Xray Grad %d w=%16.8f b=%16.8f: %16.8f %16.8f %16.8f",
//            i, width, bfactor, grad[0], grad[1], grad[2]));
//      }

      if (refinementMode.includesBFactors() || refinementMode.includesOccupancies()) {
        for (int i = nXYZ; i < n; i++) {
          g[i] = weight * gExperiment[i];
        }
      }
    } else if (refinementMode.includesBFactors() || refinementMode.includesOccupancies()) {
      // Compute the X-ray target energy and gradient.
      e = dataEnergy.energyAndGradient(x, g);
    }

    scaleCoordinatesAndGradient(x, g);
    totalEnergy = e;
    return e;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getAcceleration(double[] acceleration) {
    return dataEnergy.getAcceleration(acceleration);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getCoordinates(double[] parameters) {
    return dataEnergy.getCoordinates(parameters);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setCoordinates(double[] parameters) {
    dataEnergy.setCoordinates(parameters);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public Crystal getCrystal() {
    return dataEnergy.getCrystal();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setCrystal(Crystal crystal) {
    logger.severe(" RefinementEnergy does implement setCrystal yet.");
  }

  /**
   * Getter for the field <code>dataEnergy</code>.
   *
   * @return a {@link ffx.crystal.CrystalPotential} object.
   */
  public CrystalPotential getDataEnergy() {
    return dataEnergy;
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
    for (MolecularAssembly molecularAssembly : molecularAssemblies) {
      ForceFieldEnergy fe = molecularAssembly.getPotentialEnergy();
      fe.setEnergyTermState(state);
    }
    dataEnergy.setEnergyTermState(state);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getLambda() {
    double lambda = 1.0;
    if (data instanceof DiffractionData) {
      XRayEnergy xRayEnergy = (XRayEnergy) dataEnergy;
      lambda = xRayEnergy.getLambda();
    } else if (data instanceof RealSpaceData) {
      RealSpaceEnergy realSpaceEnergy = (RealSpaceEnergy) dataEnergy;
      lambda = realSpaceEnergy.getLambda();
    }
    return lambda;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setLambda(double lambda) {
    for (MolecularAssembly molecularAssembly : molecularAssemblies) {
      ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
      forceFieldEnergy.setLambda(lambda);
    }
    if (data instanceof DiffractionData) {
      XRayEnergy xRayEnergy = (XRayEnergy) dataEnergy;
      xRayEnergy.setLambda(lambda);
    } else if (data instanceof RealSpaceData) {
      RealSpaceEnergy realSpaceEnergy = (RealSpaceEnergy) dataEnergy;
      realSpaceEnergy.setLambda(lambda);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getMass() {
    return dataEnergy.getMass();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfVariables() {
    return dataEnergy.getNumberOfVariables();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getPreviousAcceleration(double[] previousAcceleration) {
    return dataEnergy.getPreviousAcceleration(previousAcceleration);
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
  public void setScaling(double[] scaling) {
    optimizationScaling = scaling;
  }

  /**
   * Getter for the field <code>thermostat</code>.
   *
   * @return a {@link ffx.algorithms.dynamics.thermostats.Thermostat} object.
   */
  public Thermostat getThermostat() {
    return thermostat;
  }

  /**
   * Setter for the field <code>thermostat</code>.
   *
   * @param thermostat a {@link ffx.algorithms.dynamics.thermostats.Thermostat} object.
   */
  public void setThermostat(Thermostat thermostat) {
    this.thermostat = thermostat;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalEnergy() {
    return totalEnergy;
  }

  @Override
  public List<Potential> getUnderlyingPotentials() {
    Stream<Potential> directPEs =
        Arrays.stream(molecularAssemblies).map(MolecularAssembly::getPotentialEnergy);
    Stream<Potential> allPEs =
        Arrays.stream(molecularAssemblies)
            .map(MolecularAssembly::getPotentialEnergy)
            .map(Potential::getUnderlyingPotentials)
            .flatMap(List::stream);
    return Stream.concat(directPEs, allPEs).collect(Collectors.toList());
  }

  /**
   * {@inheritDoc}
   *
   * <p>Return a reference to each variables type.
   */
  @Override
  public VARIABLE_TYPE[] getVariableTypes() {
    return dataEnergy.getVariableTypes();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getVelocity(double[] velocity) {
    return dataEnergy.getVelocity(velocity);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getd2EdL2() {
    double d2EdL2 = 0.0;
    if (thermostat != null) {
      kTScale = KCAL_TO_GRAM_ANG2_PER_PS2 / (thermostat.getTargetTemperature() * kB);
    }
    int assemblysize = molecularAssemblies.length;

    // Compute the chemical energy and gradient.
    for (int i = 0; i < assemblysize; i++) {
      ForceFieldEnergy forceFieldEnergy = molecularAssemblies[i].getPotentialEnergy();
      double curE = forceFieldEnergy.getd2EdL2();
      d2EdL2 += (curE - d2EdL2) / (i + 1);
    }
    d2EdL2 *= kTScale;

    // No 2nd derivative for scattering term.
    return d2EdL2;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getdEdL() {
    double dEdL = 0.0;
    if (thermostat != null) {
      kTScale = KCAL_TO_GRAM_ANG2_PER_PS2 / (thermostat.getTargetTemperature() * kB);
    }
    int assemblysize = molecularAssemblies.length;

    // Compute the chemical energy and gradient.
    for (int i = 0; i < assemblysize; i++) {
      ForceFieldEnergy forceFieldEnergy = molecularAssemblies[i].getPotentialEnergy();
      double curdEdL = forceFieldEnergy.getdEdL();
      dEdL += (curdEdL - dEdL) / (i + 1);
    }
    dEdL *= kTScale;
    double weight = data.getWeight();
    if (data instanceof DiffractionData) {
      XRayEnergy xRayEnergy = (XRayEnergy) dataEnergy;
      dEdL += weight * xRayEnergy.getdEdL();
    } else if (data instanceof RealSpaceData) {
      RealSpaceEnergy realSpaceEnergy = (RealSpaceEnergy) dataEnergy;
      dEdL += weight * realSpaceEnergy.getdEdL();
    }
    return dEdL;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void getdEdXdL(double[] gradient) {
    double weight = data.getWeight();
    if (thermostat != null) {
      kTScale = KCAL_TO_GRAM_ANG2_PER_PS2 / (thermostat.getTargetTemperature() * kB);
    }
    int assemblysize = molecularAssemblies.length;

    // Compute the chemical energy and gradient.
    for (int i = 0; i < assemblysize; i++) {
      ForceFieldEnergy forcefieldEnergy = molecularAssemblies[i].getPotentialEnergy();
      Arrays.fill(gChemical[i], 0.0);
      forcefieldEnergy.getdEdXdL(gChemical[i]);
    }
    for (int i = 0; i < assemblysize; i++) {
      for (int j = 0; j < nXYZ; j++) {
        gradient[j] += gChemical[i][j];
      }
    }

    // Normalize gradients for multiple-counted atoms.
    if (assemblysize > 1) {
      for (int i = 0; i < nXYZ; i++) {
        gradient[i] /= assemblysize;
      }
    }
    for (int i = 0; i < nXYZ; i++) {
      gradient[i] *= kTScale;
    }

    // Compute the X-ray target energy and gradient.
    if (gExperiment == null || gExperiment.length != nXYZ) {
      gExperiment = new double[nXYZ];
    } else {
      for (int j = 0; j < nXYZ; j++) {
        gExperiment[j] = 0.0;
      }
    }
    if (data instanceof DiffractionData) {
      XRayEnergy xRayEnergy = (XRayEnergy) dataEnergy;
      xRayEnergy.getdEdXdL(gExperiment);
    } else if (data instanceof RealSpaceData) {
      RealSpaceEnergy realSpaceEnergy = (RealSpaceEnergy) dataEnergy;
      realSpaceEnergy.getdEdXdL(gExperiment);
    }

    // Add the chemical and X-ray gradients.
    for (int i = 0; i < nXYZ; i++) {
      gradient[i] += weight * gExperiment[i];
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setAcceleration(double[] acceleration) {
    dataEnergy.setAcceleration(acceleration);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setPreviousAcceleration(double[] previousAcceleration) {
    dataEnergy.setPreviousAcceleration(previousAcceleration);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setVelocity(double[] velocity) {
    dataEnergy.setVelocity(velocity);
  }
}
