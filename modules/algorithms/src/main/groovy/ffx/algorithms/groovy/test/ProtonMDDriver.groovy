//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
import ffx.potential.bonded.AminoAcidUtils
import ffx.potential.bonded.Residue
import ffx.potential.bonded.AminoAcidUtils
import ffx.potential.cli.WriteoutOptions
import ffx.potential.extended.ExtendedSystem
import ffx.potential.extended.ExtendedVariable
import ffx.potential.extended.TautomerESV
import ffx.potential.extended.TitrationESV
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
class ProtonMDDriver extends AlgorithmsScript {

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
  ProtonMDDriver() {
    this(new Binding())
  }

  /**
   * Dynamics Constructor.
   * @param binding The Groovy Binding to use.
   */
  ProtonMDDriver(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  ProtonMDDriver run() {

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
    List<Residue> titrating = TitrationUtils.chooseTitratables(activeAssembly)
    potential = activeAssembly.getPotentialEnergy()
    double currentTemperature = dynamics.temperature
    double pe = 0
    double ke = 0
    double totalEnergy = 0

    ExtendedSystem esvSystem = new ExtendedSystem(activeAssembly)
    List<ExtendedVariable> titratingESVs = new ArrayList<>()
    List<ExtendedVariable> tautomerESVs = new ArrayList<>()
    esvSystem.setConstantPh(pH)
    for (Residue res : titrating) {
      ffx.potential.bonded.MultiResidue multi =
          TitrationUtils.titratingMultiresidueFactory(activeAssembly, res)
      TitrationESV esv = new TitrationESV(esvSystem, multi)
      titratingESVs.add(esv)
      for (Residue background : multi.getInactive()) {
        inactivateResidue(background)
      }
      esvSystem.addVariable(esv)
      if(esvSystem.config.tautomer){
        AminoAcidUtils.AminoAcid3 currentAA3 = AminoAcidUtils.AminoAcid3.valueOf(res.getName());
        if (currentAA3 == AminoAcidUtils.AminoAcid3.HIS || currentAA3 == AminoAcidUtils.AminoAcid3.HID
                || currentAA3 == AminoAcidUtils.AminoAcid3.HIE || currentAA3 == AminoAcidUtils.AminoAcid3.ASH
                || currentAA3 == AminoAcidUtils.AminoAcid3.ASP || currentAA3 == AminoAcidUtils.AminoAcid3.GLH
                || currentAA3 == AminoAcidUtils.AminoAcid3.GLU){
          TautomerESV tautomerESV = new TautomerESV(esvSystem, multi);
          tautomerESVs.add(tautomerESV);
          esvSystem.addVariable(tautomerESV);
        }
      }
    }

    ForceFieldEnergy forceFieldEnergy = (ForceFieldEnergy) potential
    forceFieldEnergy.attachExtendedSystem(esvSystem)
    logger.info(format(" Extended system with %d residues.", titratingESVs.size()))

    //Initialize thetaMass from ExtendedSystem
    int numberOfESVVariables = titratingESVs.size()
    double[] theta_mass = new double[numberOfESVVariables]
    for (int i = 0; i < numberOfESVVariables; i++) {
      theta_mass[i] = 1
    }

    // Initialize theta based on starting lambda
    double[] theta = new double[numberOfESVVariables]
    for (int i = 0; i < numberOfESVVariables; i++) {
      theta[i] = Math.asin(Math.sqrt(esvSystem.getEsv(i).getLambda()))
      logger.info(format("i: %d , lambda: %g", i, esvSystem.getEsv(i).getLambda()))
      logger.info(format("i: %d , theta: %g", i, theta[i]))
    }

    //Read in esvSystem derivatives
    double[] dEdL = new double[numberOfESVVariables]

    Random random = new Random()

    //Initialize theta velocity to Maxwell distribution based on temperature in radians/psec
    double[] theta_v = new double[numberOfESVVariables]
    for (int i = 0; i < numberOfESVVariables; i++) {
      //theta_v[i] = random.nextGaussian() * sqrt(kB * 298.15 / theta_mass[i])
      theta_v[i] = sqrt(kB * 298.15 / theta_mass[i])
    }

    double[] x = new double[potential.getNumberOfVariables()]
    double[] g = new double[potential.getNumberOfVariables()]
    potential.getCoordinates(x)
    double initialEnergy = potential.energyAndGradient(x, g, true)

    //Initialize acceleration based on gradient in radians/psec^2
    double[] esvDerivs = esvSystem.derivatives
    double[] theta_a = new double[numberOfESVVariables]
    for (int i = 0; i < numberOfESVVariables; i++) {
      dEdL[i] = esvDerivs[i] * sin(2 * theta[i])
      theta_a[i] = -KCAL_TO_GRAM_ANG2_PER_PS2 * dEdL[i] / theta_mass[i]
    }

    stochasticIntegrator = new Stochastic((double) 0.0, numberOfESVVariables, theta, theta_v,
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

      //Update ESV lambda from new theta values
      for (int i = 0; i < numberOfESVVariables; i++) {
        double sinTheta = sin(theta[i])
        //esvSystem.setLambda(i, sinTheta * sinTheta)
        esvSystem.getEsv(i).updateLambda(sinTheta * sinTheta, true)
        esvSystem.updateListeners()
      }

      //Gather derivatives
      potential.getCoordinates(x)
      pe = potential.energyAndGradient(x, g, true)
      esvDerivs = esvSystem.derivatives

      //Put derivatives in terms of theta
      for (int i = 0; i < numberOfESVVariables; i++) {
        dEdL[i] = esvDerivs[i] * sin(2 * theta[i])
      }
      //Do full step integration operation
      stochasticIntegrator.postForce(dEdL)

      //Update Kinetic Energy
      ke = 0
      for (int i = 0; i < numberOfESVVariables; i++) {
        //Compute temperature of theta particles
        double v2 = theta_v[i] * theta_v[i]
        double e = theta_mass[i] * v2
        ke += e
      }

      // Log the current state every printFrequency steps.
      currentTemperature = ke / (kB * numberOfESVVariables)
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
