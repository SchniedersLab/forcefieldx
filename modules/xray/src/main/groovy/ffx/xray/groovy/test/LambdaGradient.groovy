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
package ffx.xray.groovy.test

import org.apache.commons.configuration2.CompositeConfiguration

import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.LambdaInterface
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.GradientOptions
import ffx.xray.DiffractionData
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize.RefinementMode
import ffx.xray.cli.XrayOptions
import ffx.xray.parsers.DiffractionFile

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

/**
 * The X-ray test Lambda Gradient script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.test.LambdaGradient [options] &lt;filename&gt;
 */
@Command(description = " Test Lambda Derivatives on an X-ray target.", name = "ffxc xray.test.LambdaGradient")
class LambdaGradient extends AlgorithmsScript {

    @Mixin
    XrayOptions xrayOptions

    @Mixin
    AlchemicalOptions alchemicalOptions

    @Mixin
    GradientOptions gradientOptions

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
    private List<String> filenames
    private RefinementEnergy refinementEnergy;

    @Override
    LambdaGradient run() {

        if (!init()) {
            return this
        }

        xrayOptions.init()

        // Turn on computation of lambda derivatives
        System.setProperty("lambdaterm", "true")

        String modelfilename
        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = algorithmFunctions.openAll(filenames.get(0))
            activeAssembly = assemblies[0]
            modelfilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        } else {
            modelfilename = activeAssembly.getFile().getAbsolutePath()
            assemblies = [activeAssembly]
        }

        alchemicalOptions.setFirstSystemAlchemistry(activeAssembly)
        alchemicalOptions.setFirstSystemUnchargedAtoms(activeAssembly)

        logger.info("\n Testing X-ray lambda derivatives for " + modelfilename)

        // Load parsed X-ray properties.
        CompositeConfiguration properties = assemblies[0].getProperties()
        xrayOptions.setProperties(parseResult, properties)

        // Set up diffraction data (can be multiple files)
        List<DiffractionData> diffractionFiles = xrayOptions.processData(filenames, assemblies)
        DiffractionData diffractionData = new DiffractionData(assemblies, properties,
                xrayOptions.solventModel, diffractionFiles.toArray(new DiffractionFile[diffractionFiles.size()]))

        diffractionData.scaleBulkFit()
        diffractionData.printStats()

        refinementEnergy = new RefinementEnergy(diffractionData, xrayOptions.refinementMode)
        Potential potential = refinementEnergy
        LambdaInterface lambdaInterface = refinementEnergy

        // First atom to test.
        int atomID = gradientOptions.atomID - 1

        // Finite-difference step size in Angstroms.
        double step = gradientOptions.dx

        int n = refinementEnergy.getNumberOfVariables()
        Atom[] atoms = refinementEnergy.getActiveAtoms()
        int nAtoms = atoms.length
        double[] x = new double[n]
        double[] gradient = new double[n]

        // Finite-difference step size.
        double width = 2.0 * step
        // Error tolerence
        double errTol = 1.0e-3
        // Upper bound for typical gradient sizes (expected gradient)
        double expGrad = 1000.0

        double[] lambdaGrad = new double[n]
        double[][] lambdaGradFD = new double[2][n]

        double initialLambda = alchemicalOptions.initialLambda

        // Compute the Lambda = 0.0 energy.
        double lambda = 0.0
        lambdaInterface.setLambda(lambda)
        potential.getCoordinates(x)
        double e0 = potential.energy(x, true)

        // Compute the Lambda = 1.0 energy.
        lambda = 1.0
        lambdaInterface.setLambda(lambda)
        double e1 = potential.energy(x, true)

        logger.info(String.format(" E(0):      %20.8f.", e0))
        logger.info(String.format(" E(1):      %20.8f.", e1))
        logger.info(String.format(" E(1)-E(0): %20.8f.\n", e1 - e0))

        // Test Lambda gradient in the neighborhood of the lambda variable.
        for (int j = 0; j < 3; j++) {
            lambda = initialLambda - 0.01 + 0.01 * j

            if (lambda - step < 0.0) {
                continue
            }
            if (lambda + step > 1.0) {
                continue
            }

            logger.info(String.format(" Current lambda value %6.4f", lambda))
            lambdaInterface.setLambda(lambda)

            // Calculate the energy, dE/dX, dE/dL, d2E/dL2 and dE/dL/dX
            double e = potential.energyAndGradient(x, gradient)

            // Analytic dEdL, d2E/dL2 and dE/dL/dX
            double dEdL = lambdaInterface.getdEdL()
            double d2EdL2 = lambdaInterface.getd2EdL2()
            for (int i = 0; i < n; i++) {
                lambdaGrad[i] = 0.0
            }
            potential.getdEdXdL(lambdaGrad)

            // Calculate the finite-difference dEdLambda, d2EdLambda2 and dEdLambdadX
            lambdaInterface.setLambda(lambda + step)
            double lp = potential.energyAndGradient(x, lambdaGradFD[0])
            double dedlp = lambdaInterface.getdEdL()
            lambdaInterface.setLambda(lambda - step)
            double lm = potential.energyAndGradient(x, lambdaGradFD[1])
            double dedlm = lambdaInterface.getdEdL()

            double dEdLFD = (lp - lm) / width
            double d2EdL2FD = (dedlp - dedlm) / width

            double err = Math.abs(dEdLFD - dEdL)
            if (err < errTol) {
                logger.info(String.format(" dE/dL passed:   %10.6f", err))
            } else {
                logger.info(String.format(" dE/dL failed: %10.6f", err))
            }
            logger.info(String.format(" Numeric:   %15.8f", dEdLFD))
            logger.info(String.format(" Analytic:  %15.8f", dEdL))

            err = Math.abs(d2EdL2FD - d2EdL2)
            if (err < errTol) {
                logger.info(String.format(" d2E/dL2 passed: %10.6f", err))
            } else {
                logger.info(String.format(" d2E/dL2 failed: %10.6f", err))
            }
            logger.info(String.format(" Numeric:   %15.8f", d2EdL2FD))
            logger.info(String.format(" Analytic:  %15.8f", d2EdL2))

            boolean passed = true

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

                double error = Math.sqrt(eX * eX + eY * eY + eZ * eZ)
                if (error < errTol) {
                    logger.fine(String.format(" dE/dX/dL for Atom %d passed: %10.6f", i + 1, error))
                } else {
                    logger.info(String.format(" dE/dX/dL for Atom %d failed: %10.6f", i + 1, error))
                    logger.info(String.format(" Analytic: (%15.8f, %15.8f, %15.8f)", dXa, dYa, dZa))
                    logger.info(String.format(" Numeric:  (%15.8f, %15.8f, %15.8f)", dX, dY, dZ))
                    passed = false
                }
            }
            if (passed) {
                logger.info(String.format(" dE/dX/dL passed for all atoms"))
            }

            logger.info("")
        }

        boolean loopPrint = gradientOptions.verbose
        refinementEnergy.getCoordinates(x)
        refinementEnergy.energyAndGradient(x, gradient, loopPrint)

        double[] numeric = new double[3]
        double avLen = 0.0
        int nFailures = 0
        double avGrad = 0.0

        for (int i = atomID; i < nAtoms; i++) {
            int i3 = i * 3
            int i0 = i3 + 0
            int i1 = i3 + 1
            int i2 = i3 + 2

            // Find numeric dX
            double orig = x[i0]
            x[i0] = x[i0] + step
            double e = refinementEnergy.energy(x, loopPrint)
            x[i0] = orig - step
            e -= refinementEnergy.energy(x, loopPrint)
            x[i0] = orig
            numeric[0] = e / width

            // Find numeric dY
            orig = x[i1]
            x[i1] = x[i1] + step
            e = refinementEnergy.energy(x, loopPrint)
            x[i1] = orig - step
            e -= refinementEnergy.energy(x, loopPrint)
            x[i1] = orig
            numeric[1] = e / width

            // Find numeric dZ
            orig = x[i2]
            x[i2] = x[i2] + step
            e = refinementEnergy.energy(x, loopPrint)
            x[i2] = orig - step
            e -= refinementEnergy.energy(x, loopPrint)
            x[i2] = orig
            numeric[2] = e / width

            double dx = gradient[i0] - numeric[0]
            double dy = gradient[i1] - numeric[1]
            double dz = gradient[i2] - numeric[2]
            double len = dx * dx + dy * dy + dz * dz
            avLen += len
            len = Math.sqrt(len)

            double grad2 = gradient[i0] * gradient[i0] + gradient[i1] * gradient[i1] + gradient[i2] * gradient[i2]
            avGrad += grad2
            grad2 = Math.sqrt(grad2)

            if (len > errTol) {
                logger.info(String.format(" Atom %d failed: %10.6f.", i + 1, len)
                        + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)", gradient[i0], gradient[i1], gradient[i2])
                        + String.format("\n Numeric:  (%12.4f, %12.4f, %12.4f)", numeric[0], numeric[1], numeric[2]))
                ++nFailures
                //return
            } else {
                logger.info(String.format(" Atom %d passed: %10.6f.", i + 1, len)
                        + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)", gradient[i0], gradient[i1], gradient[i2])
                        + String.format("\n Numeric:  (%12.4f, %12.4f, %12.4f)", numeric[0], numeric[1], numeric[2]))
            }

            if (grad2 > expGrad) {
                logger.info(String.format(" Atom %d has an unusually large gradient: %10.6f", i + 1, grad2))
            }
            logger.info("\n")
        }

        avLen = avLen / nAtoms
        avLen = Math.sqrt(avLen)
        if (avLen > errTol) {
            logger.info(String.format(" Test failure: RMSD from analytic solution is %10.6f > %10.6f", avLen, errTol))
        } else {
            logger.info(String.format(" Test success: RMSD from analytic solution is %10.6f < %10.6f", avLen, errTol))
        }
        logger.info(String.format(" Number of atoms failing gradient test: %d", nFailures))

        avGrad = avGrad / nAtoms
        avGrad = Math.sqrt(avGrad)
        if (avGrad > expGrad) {
            logger.info(String.format(" Unusually large RMS gradient: %10.6f > %10.6f", avGrad, expGrad))
        } else {
            logger.info(String.format(" RMS gradient: %10.6f", avGrad))
        }

        refinementEnergy = new RefinementEnergy(diffractionData, RefinementMode.BFACTORS)
        n = refinementEnergy.getNumberOfVariables()
        gradient = new double[n]
        x = new double[n]

        refinementEnergy.getCoordinates(x)
        refinementEnergy.energyAndGradient(x, gradient)

        avLen = 0.0
        nFailures = 0
        avGrad = 0.0
        width = 2.0 * step
        errTol = 1.0e-3
        expGrad = 1000.0

        for (int i = 0; i < n; i++) {

            // Find numeric dB
            double orig = x[i]
            x[i] = x[i] + step
            double e = refinementEnergy.energy(x)
            x[i] = orig - step
            e -= refinementEnergy.energy(x)
            x[i] = orig
            double fd = e / width

            double dB = gradient[i] - fd
            double len = dB * dB
            avLen += len
            len = Math.sqrt(len)

            double grad2 = dB * dB
            avGrad += grad2
            grad2 = Math.sqrt(grad2)

            if (len > errTol) {
                logger.info(String.format(" B-Factor %d failed: %10.6f.", i + 1, len)
                        + String.format("\n Analytic: %12.4f", gradient[i])
                        + String.format("\n Numeric:  %12.4f", fd))
                ++nFailures
                //return
            } else {
                logger.info(String.format(" B-Factor %d passed: %10.6f.", i + 1, len)
                        + String.format("\n Analytic: %12.4f", gradient[i])
                        + String.format("\n Numeric:  %12.4f", fd))
            }

            if (grad2 > expGrad) {
                logger.info(String.format(" B-Factor %d has an unusually large gradient: %10.6f", i + 1, grad2))
            }
            logger.info("\n")
        }

        avLen = avLen / n
        avLen = Math.sqrt(avLen)
        if (avLen > errTol) {
            logger.info(String.format(" Test failure: RMSD from analytic solution is %10.6f > %10.6f", avLen, errTol))
        } else {
            logger.info(String.format(" Test success: RMSD from analytic solution is %10.6f < %10.6f", avLen, errTol))
        }
        logger.info(String.format(" Number of B-Factors failing gradient test: %d", nFailures))

        avGrad = avGrad / n
        avGrad = Math.sqrt(avGrad)
        if (avGrad > expGrad) {
            logger.info(String.format(" Unusually large RMS gradient: %10.6f > %10.6f", avGrad, expGrad))
        } else {
            logger.info(String.format(" RMS gradient: %10.6f", avGrad))
        }

        return this
    }

    @Override
    List<Potential> getPotentials() {
        return refinementEnergy == null ? Collections.emptyList() : Collections.singletonList(refinementEnergy);
    }
}


