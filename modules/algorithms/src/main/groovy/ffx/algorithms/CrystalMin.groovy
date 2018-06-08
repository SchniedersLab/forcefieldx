
package ffx.algorithms

import org.apache.commons.io.FilenameUtils
import static org.apache.commons.math3.util.FastMath.abs

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder

import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.MolecularAssembly.FractionalMode
import ffx.potential.XtalEnergy

/**
 * The CrystalMin script uses a limited-memory BFGS algorithm to minimize the
 * energy of a crystal, including both coordinates and unit cell parameters.
 * <br>
 * Usage:
 * <br>
 * ffxc CrystalMin [options] &lt;filename&gt;
 */
class CrystalMin extends Script {

    /**
     * Options for the Minimize Script.
     * <br>
     * Usage:
     * <br>
     * ffxc CrystalMin [options] &lt;filename&gt;
     */
    class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName='h', defaultValue='false', description='Print this help message.') boolean help;
        /**
         * -c or --coords to cycle between lattice and coordinate optimization until both satisfy the convergence criteria.
         */
        @Option(shortName='c', defaultValue='false', description='Cycle between lattice and coordinate optimization until both satisfy the convergence criteria.') boolean coords;
        /**
         * -e or --eps sets the RMS gradient convergence criteria (in kcal/mol/A^2)
         */
        @Option(shortName='e', longName='eps', defaultValue='1.0', description='RMS gradient convergence criteria') double eps;
        /**
         * -f or --fractional to set the optimization to maintain fractional coordinates [ATOM/MOLECULE/OFF].
         */
        @Option(shortName='f', defaultValue='molecule', description='Maintain fractional coordinates during lattice optimization [OFF/MOLECULE/ATOM]') String fractional;
        /**
         * -t or --tensor to print out partial derivatives of the energy with respect to unit cell parameters.
         */
        @Option(shortName='t', defaultValue='false', description='Compute partial derivatives of the energy with respect to unit cell parameters') boolean tensor;
        /**
         * The final argument should be a filename.
         */
        @Unparsed List<String> filenames;
    }

    def run() {

        def cli = new CliBuilder(usage: ' ffxc CrystalMin [options] <filename>', header: ' Options:');

        def options = new Options();
        cli.parseFromInstance(options, args);

        if (options.help == true) {
            return cli.usage();
        }

        try {
            // getAlgorithmUtils is a magic variable/closure passed in from ModelingShell
            aFuncts = getAlgorithmUtils();
        } catch (MissingMethodException ex) {
            // This is the fallback, which does everything necessary without magic names
            aFuncts = new AlgorithmUtils();
        }

        List<String> arguments = options.filenames;
        String filename = arguments.get(0);
        double eps = options.eps;

        logger.info("\n Running CrystalMinimize on " + filename);
        logger.info(" RMS gradient convergence criteria: " + eps);

        MolecularAssembly molecularAssembly = (MolecularAssembly) aFuncts.open(filename);

        ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
        XtalEnergy xtalEnergy = new XtalEnergy(forceFieldEnergy, molecularAssembly);
        xtalEnergy.setFractionalCoordinateMode(FractionalMode.MOLECULE);

        // Apply fractional coordinate mode.
        if (options.fractional) {
            try {
                FractionalMode mode = FractionalMode.valueOf(options.fractional.toUpperCase());
                xtalEnergy.setFractionalCoordinateMode(mode);
            } catch (Exception e) {
                logger.info(" Unrecognized fractional coordinate mode: " + options.fractional);
                logger.info(" Fractional coordinate mode is set to MOLECULE.");
            }
        }

        CrystalMinimize crystalMinimize = new CrystalMinimize(molecularAssembly, xtalEnergy, sh);
        crystalMinimize.minimize(eps);
        double energy = crystalMinimize.getEnergy();

        double tolerance = 1.0e-10;

        // Complete rounds of coordinate and lattice optimization.
        if (options.coords) {
            Minimize minimize = new Minimize(molecularAssembly, forceFieldEnergy, sh);
            while (true) {
                // Complete a round of coordinate optimization.
                minimize.minimize(options.eps);
                double newEnergy = minimize.getEnergy();
                double status = minimize.getStatus();
                if (status != 0) {
                    break;
                }
                if (abs(newEnergy - energy) <= tolerance) {
                    break;
                }
                energy = newEnergy;

                // Complete a round of lattice optimization.
                crystalMinimize.minimize(options.eps);
                newEnergy = crystalMinimize.getEnergy();
                status = crystalMinimize.getStatus();
                if (status != 0) {
                    break;
                }
                if (abs(newEnergy - energy) <= tolerance) {
                    break;
                }
                energy = newEnergy;
            }
        }

        if (options.tensor) {
            crystalMinimize.printTensor();
        }

        String ext = FilenameUtils.getExtension(filename);
        filename = FilenameUtils.removeExtension(filename);

        if (ext.toUpperCase().contains("XYZ")) {
            saveAsXYZ(new File(filename + ".xyz"));
        } else {
            saveAsPDB(systems, new File(filename + ".pdb"));
        }
    }
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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