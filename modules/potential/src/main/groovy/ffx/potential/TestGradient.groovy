
package ffx.potential

import groovy.cli.Option
import groovy.cli.Unparsed

import ffx.potential.bonded.Atom
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils

/**
 * The TestGradient script evaluates the consistency of the energy and gradient.
 * <br>
 * Usage:
 * <br>
 * ffxc TestGradient [options] &lt;filename&gt;
 */
class TestGradient extends Script {


    /**
     * Options for the TestGradient script.
     * <br>
     * Usage:
     * <br>
     * ffxc TestGradient [options] &lt;filename&gt;
     */
    public class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName='h', description='Print this help message.') boolean help
        /**
         * -a or --atomID to specify the number of the first atom to test.
         */
        @Option(shortName='a', defaultValue='1', description='Number of the first atom to test.') Integer atomID
        /**
         * -d or --dx to set the finite-difference step size (Angstroms).
         */
        @Option(shortName='d', defaultValue='1.0e-5', description='Finite-difference step size (Angstroms)') Double dx
        /**
         * -v or --verbose to print out the energy for each step.
         */
        @Option(shortName='v', description='Print out the energy for each step') boolean verbose
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames
    }

    /**
     * Execute the script.
     */
    def run() {

        // Create the command line parser.
        def cli = new CliBuilder(usage:' ffxc testGradient [options] <filename>', header:' Options:');
        def options = new Options()
        cli.parseFromInstance(options, args)
        if (options.help) {
            return cli.usage()
        }

        List<String> arguments = options.filenames
        String modelFilename = null
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelFilename = arguments.get(0)
        } else if (active == null) {
            return cli.usage()
        } else {
            modelFilename = active.getFile()
        }

        logger.info("\n Testing the atomic coordinate gradient of " + modelFilename + "\n");

        // This is an interface specifying the closure-like methods.
        PotentialsFunctions functions
        try {
            // Use a method closure to try to get an instance of UIUtils (the User Interfaces
            // implementation, which interfaces with the GUI, etc.).
            functions = getPotentialsUtils()
        } catch (MissingMethodException ex) {
            // If Groovy can't find the appropriate closure, catch the exception and build
            // an instance of the local implementation.
            functions = new PotentialsUtils()
        }
        // Use PotentialsFunctions methods instead of Groovy method closures to do work.
        MolecularAssembly[] assemblies = functions.open(modelFilename)
        MolecularAssembly activeAssembly = assemblies[0]
        ForceFieldEnergy energy = activeAssembly.getPotentialEnergy();
        Atom[] atoms = activeAssembly.getAtomArray();
        int n = atoms.length;

        // First atom to test.
        int atomID = options.atomID - 1
        if (atomID >= n) {
            atomID = 0;
        }
        logger.info("\n First atom to test:\t\t" + (atomID + 1))

        // Finite-difference step size in Angstroms.
        double step = options.dx
        logger.info(" Finite-difference step size:\t" + step);

        // Print out the energy for each step.
        boolean print = options.verbose
        logger.info(" Verbose printing:\t\t" + print + "\n");

        // Upper bound for typical gradient sizes (expected gradient)
        double expGrad = 1000.0;
        double gradientTolerance = 1.0e-3;
        double width = 2.0 * step;
        double[] x = new double[n*3];
        double[] analytic = new double[3];
        double[] numeric = new double[3];
		
        double avLen = 0.0;
        int nFailures = 0;
        double avGrad = 0.0;
        double expGrad2 = expGrad * expGrad;
        for (int i=atomID; i<n; i++) {
            Atom a0 = atoms[i];
            int i3 = i*3;
            int i0 = i3 + 0;
            int i1 = i3 + 1;
            int i2 = i3 + 2;
			energy.energy(true, true);
			a0.getXYZGradient(analytic);
			energy.getCoordinates(x);

            // Find numeric dX
            double orig = x[i0];
            x[i0] = x[i0] + step;
            energy.setCoordinates(x);
            double e = energy.energy(false, print);
            x[i0] = orig - step;
            energy.setCoordinates(x);
            e -= energy.energy(false, print);
            x[i0] = orig;
            numeric[0] = e / width;

            // Find numeric dY
            orig = x[i1];
            x[i1] = x[i1] + step;
            energy.setCoordinates(x);
            e = energy.energy(false, print);
            x[i1] = orig - step;
            energy.setCoordinates(x);
            e -= energy.energy(false, print);
            x[i1] = orig;
            numeric[1] = e / width;

            // Find numeric dZ
            orig = x[i2];
            x[i2] = x[i2] + step;
            energy.setCoordinates(x);
            e = energy.energy(false, print);
            x[i2] = orig - step;
            energy.setCoordinates(x);
            e -= energy.energy(false, print);
            x[i2] = orig;
            numeric[2] = e / width;

            double dx = analytic[0] - numeric[0];
            double dy = analytic[1] - numeric[1];
            double dz = analytic[2] - numeric[2];
            double len = dx * dx + dy * dy + dz * dz;
            avLen += len;
            len = Math.sqrt(len);

            double grad2 = analytic[0] * analytic[0] + analytic[1] * analytic[1] + analytic[2] * analytic[2];
            avGrad += grad2;

            if (len > gradientTolerance) {
                logger.info(" " + a0.toShortString() + String.format(" failed: %10.6f.", len)
                    + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[0], analytic[1], analytic[2])
                    + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]));
                ++nFailures;
                //return;
            } else {
                logger.info(" " + a0.toShortString() + String.format(" passed: %10.6f.", len)
                    + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[0], analytic[1], analytic[2])
                    + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)", numeric[0], numeric[1], numeric[2]));
            }

            if (grad2 > expGrad2) {
                logger.info(String.format(" Atom %d has an unusually large gradient: %10.6f", i+1, Math.sqrt(grad2)));
            }
            logger.info("\n");
        }

        avLen = avLen / n;
        avLen = Math.sqrt(avLen);
        if (avLen > gradientTolerance) {
            logger.info(String.format(" Test failure: RMSD from analytic solution is %10.6f > %10.6f", avLen, gradientTolerance));
        } else {
            logger.info(String.format(" Test success: RMSD from analytic solution is %10.6f < %10.6f", avLen, gradientTolerance));
        }
        logger.info(String.format(" Number of atoms failing analytic test: %d", nFailures));

        avGrad = avGrad / n;
        avGrad = Math.sqrt(avGrad);
        if (avGrad > expGrad) {
            logger.info(String.format(" Unusually large RMS gradient: %10.6f > %10.6f", avGrad, expGrad));
        } else {
            logger.info(String.format(" RMS gradient: %10.6f", avGrad));
        }

    }
}
/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */