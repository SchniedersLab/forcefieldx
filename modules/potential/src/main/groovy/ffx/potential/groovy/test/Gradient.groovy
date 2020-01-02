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

import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.GradientOptions
import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

/**
 * The Gradient script evaluates the consistency of the energy and gradient.
 * <br>
 * Usage:
 * <br>
 * ffxc test.Gradient [options] &lt;filename&gt;
 */
@Command(description = " Test the potential energy gradient.", name = "ffxc test.Gradient")
class Gradient extends PotentialScript {

    @Mixin
    GradientOptions gradientOptions

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    private ForceFieldEnergy energy;
    public int nFailures = 0

    /**
     * Execute the script.
     */
    @Override
    Gradient run() {

        if (!init()) {
            return this
        }

        String modelFilename
        if (filenames != null && filenames.size() > 0) {
            modelFilename = filenames.get(0)
            MolecularAssembly[] assemblies = potentialFunctions.openAll(modelFilename)
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        } else {
            modelFilename = activeAssembly.getFile().getAbsolutePath()
        }

        logger.info("\n Testing the atomic coordinate gradient of " + modelFilename + "\n");

        energy = activeAssembly.getPotentialEnergy()
        Atom[] atoms = activeAssembly.getAtomArray()
        int nAtoms = atoms.length

        // First atom to test.
        int atomID = gradientOptions.atomID - 1
        if (atomID >= nAtoms) {
            atomID = 0
        }
        logger.info("\n First atom to test:\t\t" + (atomID + 1))

        // Last atom to test.
        int lastAtomID = gradientOptions.lastAtomID - 1
        if (lastAtomID < atomID || lastAtomID >= nAtoms) {
            lastAtomID = nAtoms - 1
        }
        logger.info("\n Last atom to test:\t\t" + (lastAtomID + 1))

        // Finite-difference step size in Angstroms.
        double step = gradientOptions.dx
        logger.info(" Finite-difference step size:\t" + step)

        // Print out the energy for each step.
        boolean print = gradientOptions.verbose
        logger.info(" Verbose printing:\t\t" + print + "\n")

        // Upper bound for typical gradient sizes (expected gradient)
        double expGrad = 1000.0
        double gradientTolerance = gradientOptions.tolerance
        double width = 2.0 * step
        double[] x = new double[nAtoms * 3]
        double[] analytic = new double[3]
        double[] numeric = new double[3]

        double avLen = 0.0
        double avGrad = 0.0
        double expGrad2 = expGrad * expGrad
        for (int i = atomID; i <= lastAtomID; i++) {
            Atom a0 = atoms[i]
            int i3 = i * 3
            int i0 = i3 + 0
            int i1 = i3 + 1
            int i2 = i3 + 2
            energy.energy(true, true)
            a0.getXYZGradient(analytic)
            energy.getCoordinates(x)

            // Find numeric dX
            double orig = x[i0]
            x[i0] = x[i0] + step
            energy.setCoordinates(x)
            double e = energy.energy(false, print)
            x[i0] = orig - step
            energy.setCoordinates(x)
            e -= energy.energy(false, print)
            x[i0] = orig
            numeric[0] = e / width

            // Find numeric dY
            orig = x[i1]
            x[i1] = x[i1] + step
            energy.setCoordinates(x)
            e = energy.energy(false, print)
            x[i1] = orig - step
            energy.setCoordinates(x)
            e -= energy.energy(false, print)
            x[i1] = orig
            numeric[1] = e / width

            // Find numeric dZ
            orig = x[i2]
            x[i2] = x[i2] + step
            energy.setCoordinates(x)
            e = energy.energy(false, print)
            x[i2] = orig - step
            energy.setCoordinates(x)
            e -= energy.energy(false, print)
            x[i2] = orig
            numeric[2] = e / width

            double dx = analytic[0] - numeric[0]
            double dy = analytic[1] - numeric[1]
            double dz = analytic[2] - numeric[2]
            double len = dx * dx + dy * dy + dz * dz
            avLen += len
            len = Math.sqrt(len)

            double grad2 = analytic[0] * analytic[0] + analytic[1] * analytic[1] + analytic[2] * analytic[2]
            avGrad += grad2

            if (len > gradientTolerance) {
                logger.info(" " + a0.toShortString() + String.format(" failed: %10.6f.", len)
                        + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[0], analytic[1], analytic[2])
                        + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]))
                ++nFailures
                //return
            } else {
                logger.info(" " + a0.toShortString() + String.format(" passed: %10.6f.", len)
                        + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[0], analytic[1], analytic[2])
                        + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)", numeric[0], numeric[1], numeric[2]))
            }

            if (grad2 > expGrad2) {
                logger.info(String.format(" Atom %d has an unusually large gradient: %10.6f", i + 1, Math.sqrt(grad2)))
            }
            logger.info("\n")
        }

        avLen = avLen / nAtoms
        avLen = Math.sqrt(avLen)
        if (avLen > gradientTolerance) {
            logger.info(String.format(" Test failure: RMSD from analytic solution is %10.6f > %10.6f", avLen, gradientTolerance))
        } else {
            logger.info(String.format(" Test success: RMSD from analytic solution is %10.6f < %10.6f", avLen, gradientTolerance))
        }
        logger.info(String.format(" Number of atoms failing analytic test: %d", nFailures))

        avGrad = avGrad / nAtoms
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
        return energy == null ? Collections.emptyList() : Collections.singletonList(energy);
    }
}
