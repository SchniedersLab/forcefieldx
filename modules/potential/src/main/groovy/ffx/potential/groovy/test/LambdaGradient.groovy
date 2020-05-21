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
package ffx.potential.groovy.test

import edu.rit.pj.ParallelTeam
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergyOpenMM
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.LambdaInterface
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.GradientOptions
import ffx.potential.cli.PotentialScript
import ffx.potential.cli.TopologyOptions
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.stream.IntStream

import static ffx.utilities.StringUtils.parseAtomRanges
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.abs
import static org.apache.commons.math3.util.FastMath.sqrt

/**
 * The LambdaGradient script tests numeric gradients w.r.t. lambda against
 * numeric, finite-difference gradients
 * <br>
 * Usage:
 * <br>
 * ffxc test.LambdaGradient [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Test potential energy derivatives with respect to Lambda.", name = "ffxc test.LambdaGradient")
class LambdaGradient extends PotentialScript {

  @Mixin
  private AlchemicalOptions alchemical

  @Mixin
  private TopologyOptions topology

  @Mixin
  private GradientOptions gradientOptions

  /**
   * --ls or --lambdaScan Scan lambda values.
   */
  @Option(names = ['--ls', '--lambdaScan'], paramLabel = 'false',
      description = 'Scan lambda values.')
  boolean lambdaScan = false

  /**
   * --lm or --lambdaMoveSize Size of the lambda moves during the test.
   */
  @Option(names = ['--lm', '--lambdaMoveSize'], paramLabel = '0.01',
      description = 'Size of the lambda moves during the test.')
  double lambdaMoveSize = 0.01

  /**
   * --sk2 or --skip2 disables dU2/dL2, and dU2/dLdX tests.
   */
  @Option(names = ['--sk2', '--skip2'], paramLabel = 'false',
      description = 'Skip 2nd derivatives.')
  boolean skipSecondDerivatives = false

  /**
   * -sdX or --skipdX disables the very long atomic gradient tests and
   * only performs dU/dL, dU2/dL2, and dU2/dLdX tests. Useful if the
   * script is being used to generate a dU/dL profile.
   */
  @Option(names = ['--sdX', '--skipdX'], paramLabel = 'false',
      description = 'Skip calculating per-atom dUdX values and only test lambda gradients.')
  boolean skipAtomGradients = false

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = 'The atomic coordinate file in PDB or XYZ format.')
  List<String> filenames = null

  private int threadsAvail = ParallelTeam.getDefaultThreadCount()
  private int threadsPer = threadsAvail
  MolecularAssembly[] topologies
  private Potential potential

  public int ndEdLFailures = 0
  public int ndEdXFailures = 0
  public int nd2EdL2Failures = 0
  public int ndEdXdLFailures = 0
  public double e0 = 0.0
  public double e1 = 0.0

  /**
   * LambdaGradient Constructor
   */
  LambdaGradient() {
    this(new Binding())
  }

  /**
   * LambdaGradient Constructor
   * @param binding Groovy Binding to use.
   */
  LambdaGradient(Binding binding) {
    super(binding)
  }

  /**
   * Script run method.
   */
  @Override
  LambdaGradient run() {

    if (!init()) {
      return this
    }

    List<String> arguments = filenames
    // Check nArgs should either be number of arguments (min 1), else 1.
    int nArgs = arguments ? arguments.size() : 1
    nArgs = (nArgs < 1) ? 1 : nArgs

    topologies = new MolecularAssembly[nArgs]

    int numParallel = topology.getNumParallel(threadsAvail, nArgs)
    threadsPer = (int) (threadsAvail / numParallel)

    // Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
    /* Checking nArgs == 1 should only be done for scripts that imply some sort of lambda scaling.
       The Minimize script, for example, may be running on a single, unscaled physical topology. */
    boolean lambdaTerm = (nArgs == 1 || alchemical.hasSoftcore() || topology.hasSoftcore())

    if (lambdaTerm) {
      System.setProperty("lambdaterm", "true")
    }

    alchemical.getInitialLambda()

    // Relative free energies via the DualTopologyEnergy class require different
    // default OST parameters than absolute free energies.
    if (nArgs >= 2) {
      // Ligand vapor electrostatics are not calculated. This cancels when the
      // difference between protein and water environments is considered.
      System.setProperty("ligand-vapor-elec", "false")
    }

    List<MolecularAssembly> topologyList = new ArrayList<>(4)

    // Read in files.
    if (!arguments || arguments.isEmpty()) {
      MolecularAssembly mola = potentialFunctions.getActiveAssembly()
      if (mola == null) {
        logger.info(helpString())
        return this
      }
      arguments = new ArrayList<>()
      arguments.add(mola.getFile().getName())
      topologyList.add(alchemical.processFile(topology, mola, 0))
    } else {
      logger.info(format(" Initializing %d topologies...", nArgs))
      for (int i = 0; i < nArgs; i++) {
        topologyList.add(alchemical.openFile(potentialFunctions, topology, threadsPer,
            arguments.get(i), i))
      }
    }

    MolecularAssembly[] topologies =
        topologyList.toArray(new MolecularAssembly[topologyList.size()])

    // Configure the potential to test.
    StringBuilder sb = new StringBuilder("\n Testing lambda derivatives for ")
    potential = topology.assemblePotential(topologies, threadsAvail, sb)

    if (potential instanceof ForceFieldEnergyOpenMM) {
      skipSecondDerivatives = true
    }

    logger.info(sb.toString())

    LambdaInterface lambdaInterface = (LambdaInterface) potential

    // Reset the number of variables for the case of dual topology.
    int n = potential.getNumberOfVariables()
    double[] x = new double[n]
    double[] gradient = new double[n]
    double[] lambdaGrad = new double[n]
    double[][] lambdaGradFD = new double[2][n]

    // Number of active atoms.
    assert (n % 3 == 0)
    int nAtoms = (int) (n / 3)

    // Compute the Lambda = 0.0 energy.
    double lambda = 0.0
    lambdaInterface.setLambda(lambda)
    potential.getCoordinates(x)

    e0 = potential.energyAndGradient(x, gradient)
    double dEdL = lambdaInterface.getdEdL()
    logger.info(format(" L=%4.2f E=%12.6f dE/dL=%12.6f", lambda, e0, dEdL))

    // Scan intermediate lambda values.
    if (lambdaScan) {
      for (int i = 1; i <= 9; i++) {
        lambda = i * 0.1
        lambdaInterface.setLambda(lambda)
        double e = potential.energyAndGradient(x, gradient)
        dEdL = lambdaInterface.getdEdL()
        logger.info(format(" L=%4.2f E=%12.6f dE/dL=%12.6f", lambda, e, dEdL))
      }
    }

    // Compute the Lambda = 1.0 energy.
    lambda = 1.0
    lambdaInterface.setLambda(lambda)
    e1 = potential.energyAndGradient(x, gradient)
    dEdL = lambdaInterface.getdEdL()
    logger.info(format(" L=%4.2f E=%12.6f dE/dL=%12.6f", lambda, e1, dEdL))
    logger.info(format(" E(1)-E(0): %12.6f.\n", e1 - e0))

    // Finite-difference step size.
    double step = gradientOptions.dx
    double width = 2.0 * step
    logger.info(" Finite-difference step size:\t" + step)

    // Print out the energy for each step.
    boolean print = gradientOptions.verbose
    logger.info(" Verbose printing:\t\t" + print + "\n")

    // Error tolerence
    double errTol = gradientOptions.tolerance

    // Upper bound for typical gradient sizes (expected gradient)
    double expGrad = 1000.0

    // Test Lambda gradient in the neighborhood of the lambda variable.
    for (int j = 0; j < 3; j++) {
      // Loop-local counter for printout.
      int jd2EdXdLFailures = 0

      lambda = alchemical.initialLambda - lambdaMoveSize + lambdaMoveSize * j

      if (lambda - gradientOptions.dx < 0.0 || lambda + gradientOptions.dx > 1.0) {
        continue
      }

      logger.info(format(" Current lambda value %6.4f", lambda))
      lambdaInterface.setLambda(lambda)

      // Calculate the energy, dE/dX, dE/dL, d2E/dL2 and dE/dL/dX
      double e = potential.energyAndGradient(x, gradient)

      // Analytic dEdL, d2E/dL2 and dE/dL/dX
      dEdL = lambdaInterface.getdEdL()

      double d2EdL2 = lambdaInterface.getd2EdL2()
      for (int i = 0; i < n; i++) {
        lambdaGrad[i] = 0.0
      }
      lambdaInterface.getdEdXdL(lambdaGrad)

      // Calculate the finite-difference dEdLambda, d2EdLambda2 and dEdLambdadX
      lambdaInterface.setLambda(lambda + gradientOptions.dx)
      double lp = potential.energyAndGradient(x, lambdaGradFD[0])
      double dedlp = lambdaInterface.getdEdL()
      lambdaInterface.setLambda(lambda - gradientOptions.dx)
      double lm = potential.energyAndGradient(x, lambdaGradFD[1])
      double dedlm = lambdaInterface.getdEdL()

      double dEdLFD = (lp - lm) / width
      double d2EdL2FD = (dedlp - dedlm) / width

      double err = abs(dEdLFD - dEdL)
      if (err < errTol) {
        logger.info(format(" dE/dL passed:   %10.6f", err))
      } else {
        logger.info(format(" dE/dL failed: %10.6f", err))
        ndEdLFailures++
      }
      logger.info(format(" Numeric:   %15.8f", dEdLFD))
      logger.info(format(" Analytic:  %15.8f", dEdL))

      if (!skipSecondDerivatives) {
        err = abs(d2EdL2FD - d2EdL2)
        if (err < errTol) {
          logger.info(format(" d2E/dL2 passed: %10.6f", err))
        } else {
          logger.info(format(" d2E/dL2 failed: %10.6f", err))
          nd2EdL2Failures++
        }
        logger.info(format(" Numeric:   %15.8f", d2EdL2FD))
        logger.info(format(" Analytic:  %15.8f", d2EdL2))

        double rmsError = 0
        for (int i = 0; i < nAtoms; i++) {
          int ii = i * 3
          double dX = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / width
          double dXa = lambdaGrad[ii]
          double eX = dX - dXa
          ii++
          double dY = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / width
          double dYa = lambdaGrad[ii]
          double eY = dY - dYa
          ii++
          double dZ = (lambdaGradFD[0][ii] - lambdaGradFD[1][ii]) / width
          double dZa = lambdaGrad[ii]
          double eZ = dZ - dZa

          double error = eX * eX + eY * eY + eZ * eZ
          rmsError += error
          error = sqrt(error)
          if (error < errTol) {
            logger.fine(format(" dE/dX/dL for degree of freedom %d passed: %10.6f", i + 1, error))
          } else {
            logger.info(format(" dE/dX/dL for degree of freedom %d failed: %10.6f", i + 1, error))
            logger.info(format(" Analytic: (%15.8f, %15.8f, %15.8f)", dXa, dYa, dZa))
            logger.info(format(" Numeric:  (%15.8f, %15.8f, %15.8f)", dX, dY, dZ))
            ndEdXdLFailures++
            jd2EdXdLFailures++
          }
        }
        rmsError = sqrt(rmsError / nAtoms)
        if (ndEdXdLFailures == 0) {
          logger.info(
              format(" dE/dX/dL passed for all degrees of freedom: RMS error %15.8f", rmsError))
        } else {
          logger.info(
              format(" dE/dX/dL failed for %d of %d atoms: RMS error %15.8f", jd2EdXdLFailures,
                  nAtoms, rmsError))
        }
        logger.info("")
      }
    }

    lambdaInterface.setLambda(alchemical.initialLambda)
    potential.getCoordinates(x)
    potential.energyAndGradient(x, gradient, print)

    if (!skipAtomGradients) {
      double[] numeric = new double[3]
      double avLen = 0.0
      double avGrad = 0.0

      // Collect atoms to test.
      List<Integer> degreesOfFreedomToTest
      if (gradientOptions.gradientAtoms.equalsIgnoreCase("NONE")) {
        logger.info(" The gradient of no atoms will be evaluated.")
        return this
      } else if (gradientOptions.gradientAtoms.equalsIgnoreCase("ALL")) {
        logger.info(" Checking gradient for all active atoms.\n")
        degreesOfFreedomToTest = new ArrayList<>()
        IntStream.range(0, nAtoms).forEach(val -> degreesOfFreedomToTest.add(val))
      } else {
        degreesOfFreedomToTest =
            parseAtomRanges(" Gradient atoms", gradientOptions.gradientAtoms, nAtoms)
        logger.info(
            " Checking gradient for active atoms in the range: " + gradientOptions.gradientAtoms
                +
                "\n")
      }

      for (int i : degreesOfFreedomToTest) {
        int i3 = i * 3
        int i0 = i3 + 0
        int i1 = i3 + 1
        int i2 = i3 + 2

        // Find numeric dX
        double orig = x[i0]
        x[i0] = x[i0] + step
        double e = potential.energyAndGradient(x, lambdaGradFD[0], print)
        x[i0] = orig - step
        e -= potential.energyAndGradient(x, lambdaGradFD[1], print)
        x[i0] = orig
        numeric[0] = e / width

        // Find numeric dY
        orig = x[i1]
        x[i1] = x[i1] + step
        e = potential.energyAndGradient(x, lambdaGradFD[0], print)
        x[i1] = orig - step
        e -= potential.energyAndGradient(x, lambdaGradFD[1], print)
        x[i1] = orig
        numeric[1] = e / width

        // Find numeric dZ
        orig = x[i2]
        x[i2] = x[i2] + step
        e = potential.energyAndGradient(x, lambdaGradFD[0], print)
        x[i2] = orig - step
        e -= potential.energyAndGradient(x, lambdaGradFD[1], print)
        x[i2] = orig
        numeric[2] = e / width

        double dx = gradient[i0] - numeric[0]
        double dy = gradient[i1] - numeric[1]
        double dz = gradient[i2] - numeric[2]
        double len = dx * dx + dy * dy + dz * dz
        avLen += len
        len = Math.sqrt(len)

        double grad2 =
            gradient[i0] * gradient[i0] + gradient[i1] * gradient[i1] + gradient[i2] * gradient[i2]
        avGrad += grad2

        if (len > errTol) {
          logger.info(format(" Degree of freedom %d failed: %10.6f.", i + 1, len)
              + format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", gradient[i0], gradient[i1],
              gradient[i2])
              + format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]))
          ++ndEdXFailures
        } else {
          logger.info(format(" Degree of freedom %d passed: %10.6f.", i + 1, len)
              + format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", gradient[i0], gradient[i1],
              gradient[i2])
              + format(" Numeric:  (%12.4f, %12.4f, %12.4f)", numeric[0], numeric[1], numeric[2]))
        }

        if (grad2 > expGrad) {
          logger.info(
              format(" Degree of freedom %d has an unusually large gradient: %10.6f", i + 1, grad2))
        }
        logger.info("\n")
      }

      avLen = avLen / nAtoms
      avLen = sqrt(avLen)
      if (avLen > errTol) {
        logger.info(
            format(" Test failure: RMSD from analytic solution is %10.6f > %10.6f", avLen, errTol))
      } else {
        logger.info(
            format(" Test success: RMSD from analytic solution is %10.6f < %10.6f", avLen, errTol))
      }
      logger.info(format(" Number of atoms failing gradient test: %d", ndEdXFailures))

      avGrad = avGrad / nAtoms
      avGrad = sqrt(avGrad)
      if (avGrad > expGrad) {
        logger.info(format(" Unusually large RMS gradient: %10.6f > %10.6f", avGrad, expGrad))
      } else {
        logger.info(format(" RMS gradient: %10.6f", avGrad))
      }
    } else {
      logger.info(" Skipping atomic dU/dX gradients.")
    }

    return this
  }

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
