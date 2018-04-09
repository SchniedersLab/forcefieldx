
import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed

import ffx.algorithms.CrystalMinimize
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.XtalEnergy


/**
 * The CrystalMinimize script uses a limited-memory BFGS algorithm to minimize the
 * energy of a crystal, including both coordinates and unit cell parameters.
 * <br>
 * Usage:
 * <br>
 * ffxc CrystalMinimize [options] &lt;filename&gt;
 */
class CrystalMinimize extends Script {

    /**
     * Options for the Minimize Script.
     * <br>
     * Usage:
     * <br>
     * ffxc CrystalMinimize [options] &lt;filename&gt;
     */
    class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName='h', defaultValue='false', description='Print this help message.') boolean help;
        /**
         * -e or --eps sets the RMS gradient convergence criteria (in kcal/mol/A^2)
         */
        @Option(shortName='e', longName='eps', defaultValue='1.0', description='RMS gradient convergence criteria') double eps;
        /**
         * The final argument should be a filename.
         */
        @Unparsed List<String> filenames;
    }

    def run() {

        def cli = new CliBuilder(usage: ' ffxc CrystalMinimize [options] <filename>', header: ' Options:');

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
            aFuncts = new ffx.algorithms.AlgorithmUtils();
        }

        List<String> arguments = options.filenames;
        String filename = arguments.get(0);
        double eps = options.eps;

        logger.info("\n Running CrystalMinimize on " + filename);
        logger.info(" RMS gradient convergence criteria: " + eps);

        MolecularAssembly molecularAssembly = aFuncts.open(filename);

        ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
        XtalEnergy xtalEnergy = new XtalEnergy(forceFieldEnergy, molecularAssembly);

        // Do the minimization
        CrystalMinimize crystalMinimize = new CrystalMinimize(molecularAssembly, xtalEnergy, sh);

        double e = crystalMinimize.minimize(eps);

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