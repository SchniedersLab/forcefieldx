//******************************************************************************
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
//******************************************************************************
package ffx.algorithms.groovy.test

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.dynamics.integrators.Stochastic
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Residue
import ffx.potential.cli.WriteoutOptions
import ffx.potential.extended.NewExtendedSystem
import ffx.potential.extended.TitrationUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static ffx.potential.extended.TitrationUtils.inactivateResidue
import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2
import static ffx.utilities.Constants.kB
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.sin
import static org.apache.commons.math3.util.FastMath.sqrt

/**
 * The Dynamics script implements molecular and stochastic dynamics algorithms.
 * <br>
 * Usage:
 * <br>
 * ffxc CPHMDynamics [options] &lt;filename&gt; [file2...]
 */
@Command(description = " Run constant pH dynamics on a system.", name = "ffxc PHDynamics")
class NewProtonMDDriver extends AlgorithmsScript {

  @Mixin
  DynamicsOptions dynamics

  @Mixin
  WriteoutOptions writeOut

  /**
   * --pH or --constantPH Constant pH value for molecular dynamics.
   */
  @Option(names = ['--pH', '--constantPH'], paramLabel = '7.4',
          description = 'Constant pH value for molecular dynamics')
  double pH = 7.4

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
          description = "XYZ or PDB input files.")
  private List<String> filenames

  /**
   * Creation of a public field to try and make the JUnit test work, original code does not declare this as a public field.
   * Originally it is declared in the run method
   */
  public Potential potential
  public Stochastic stochasticIntegrator = null
  private double totalSimTime = 0.0

  Potential getPotentialObject() {
    return potential
  }

  /**
   * Dynamics Constructor.
   */
  NewProtonMDDriver() {
    this(new Binding())
  }

  /**
   * Dynamics Constructor.
   * @param binding The Groovy Binding to use.
   */
  NewProtonMDDriver(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  NewProtonMDDriver run() {

    if (!init()) {
      return this
    }

    TitrationUtils.initDiscountPreloadProperties()

    String modelFilename
    if (filenames != null && filenames.size() > 0) {
      MolecularAssembly molecularAssembly =
              TitrationUtils.openFullyProtonated(new File(filenames.get(0)))
      MolecularAssembly[] assemblies = [molecularAssembly]
      activeAssembly = assemblies[0]
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    } else {
      modelFilename = activeAssembly.getFile().getAbsolutePath()
    }

    // Select all possible titrating residues.
    potential = activeAssembly.getPotentialEnergy()
    double currentTemperature = dynamics.temperature
    double pe = 0
    double ke = 0
    double totalEnergy = 0

    NewExtendedSystem esvSystem = new NewExtendedSystem(activeAssembly)
    esvSystem.setConstantPh(pH)
    List<Residue> extendedResidues = esvSystem.getExtendedResidueList()
    List<Residue> titratingResidues = esvSystem.getTitratingResidueList()
    List<Residue> tautomerResidues = esvSystem.getTautomerizingResidueList()

    //energy.attachExtendedSystem(esvSystem)
    int numESVs = extendedResidues.size()
    logger.info(format(" Extended system with %d residues.", numESVs))

    ForceFieldEnergy forceFieldEnergy = (ForceFieldEnergy) potential
    logger.info(format(" Extended system with %d residues.", numESVs))

    //Initialize thetaMass from ExtendedSystem
    double[] theta_mass = esvSystem.getThetaMassArray()

    // Initialize theta based on starting lambda
    double[] theta = esvSystem.getThetaPosition()
    for (int i=0; i < theta.size(); i++) {
      logger.info(format("i: %d , lambda: %g", i, esvSystem.getExtendedLambdas()[i]))
      logger.info(format("i: %d , theta: %g", i, theta[i]))
    }

    //Read in esvSystem derivatives
    double[] dEdL = new double[numESVs]

    //Initialize theta velocity to Maxwell distribution based on temperature in radians/psec
    double[] theta_v = esvSystem.getThetaVelocity()
    for (int i = 0; i < theta_v.size(); i++) {
      logger.info(format("i: %d , theta_v: %g", i, theta_v[i]))
    }

    //Initialize acceleration based on gradient in radians/psec^2
    //double[] esvDerivs = esvSystem.derivatives
    double[] theta_a = esvSystem.getThetaAccel()
    for (int i = 0; i < theta_a.size(); i++) {
      logger.info(format("i: %d , theta_a: %g", i, theta_a[i]))
    }

    double[] x = new double[potential.getNumberOfVariables()]
    double[] g = new double[potential.getNumberOfVariables()]
    potential.getCoordinates(x)
    double initialEnergy = potential.energyAndGradient(x, g, true)

    stochasticIntegrator = new Stochastic((double) 0.0, numESVs, theta, theta_v,
            theta_a, theta_mass)
    logger.info("\n Running lambda dynamics on " + modelFilename)
    logger.info(format("%s", esvSystem.getLambdaList()))

    for (long step = 1; step <= dynamics.steps; step++) {
      // Update extended system variables if present.
      /*if (esvSystem != null) {
        esvSystem.propagateESVs(dynamics.temperature, dynamics.dt, step * dynamics.dt)
      }*/
      totalSimTime += dynamics.dt

      //Do half step integration operation
      stochasticIntegrator.preForce(potential)
      esvSystem.preForce()

      //Gather derivatives
      potential.getCoordinates(x)
      pe = potential.energyAndGradient(x, g, true)
      pe += esvSystem.getBiasEnergy()
      //Put derivatives in terms of theta
      dEdL = esvSystem.postForce()
      //Do full step integration operation
      stochasticIntegrator.postForce(dEdL)

      //Update Kinetic Energy
      ke = 0
      for (int i = 0; i < numESVs; i++) {
        //Compute temperature of theta particles
        double v2 = theta_v[i] * theta_v[i]
        double e = theta_mass[i] * v2
        ke += e
      }

      // Log the current state every printFrequency steps.
      currentTemperature = ke / (kB * numESVs)
      //Degrees of freedom for theta particle on unit circle??
      ke *= 0.5 / KCAL_TO_GRAM_ANG2_PER_PS2
      totalEnergy = pe + ke
      logger.info(
              format("Potential Energy: %16.8f Kinetic Energy %16.8f Total Energy: %16.8f", pe, ke,
                      totalEnergy))
      logger.info(format("Current Temperature: %g", currentTemperature))
      logger.info(format(" %7.3e %s", totalSimTime, esvSystem.getLambdaList()))
    }

    return this
  }

  /**
   * {@inheritDoc}
   */
  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (potential == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList(potential)
    }
    return potentials
  }

}
