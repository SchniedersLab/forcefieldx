
package ffx.potential

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder;

import ffx.potential.bonded.Atom
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils

/**
 * The Energy script evaluates the energy of a system.
 * <br>
 * Usage:
 * <br>
 * ffxc Energy &lt;filename&gt;
 */
class Energy extends Script {

    /**
     * Options for the Energy Script.
     * <br>
     * Usage:
     * <br>
     * ffxc Energy [options] &lt;filename&gt;
     */
    public class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName='h', defaultValue='false', description='Print this help message.') boolean help;
        /**
         * -g or --gradient to print out gradients.
         */
        @Option(shortName='g', longName='gradient', description='Print out atomic gradients as well as energy.') boolean gradient;
        /**
         * -es1 or --noElecStart1 defines the first atom of the first topology to have no electrostatics.
         */
        @Option(shortName='es1', longName='noElecStart1', defaultValue='1', description='Starting no-electrostatics atom for 1st topology') int es1;
        /**
         * -ef1 or --noElecFinal1 defines the last atom of the first topology to have no electrostatics.
         */
        @Option(shortName='ef1', longName='noElecFinal1', defaultValue='-1', description='Final no-electrostatics atom for 1st topology') int ef1;
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed(description='The atomic coordinate file in PDB or XYZ format.') List<String> filenames
    }


    /**
     * Execute the script.
     */
    def run() {

        // def cli = new CliBuilder(usage:' ffxc Energy [options] <filename>', header:' Options:');
        def cli = new CliBuilder();
        cli.name = "ffxc Energy";

        // String args[] = getArgs();

        def options = new Options()
        cli.parseFromInstance(options, args)

        if (options.help == true) {
            return cli.usage()
        }

        List<String> arguments = options.filenames
        String modelFilename
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelFilename = arguments.get(0)
            //open(modelFilename)
        } else if (active == null) {
            return cli.usage()
        } else {
            modelFilename = active.getFile()
        }

        logger.info("\n Running Energy on " + modelFilename)

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
        MolecularAssembly activeAssembly = assemblies[0];
        ForceFieldEnergy pe = activeAssembly.getPotentialEnergy();

        Atom[] atoms = activeAssembly.getAtomArray();

        // Apply the no electrostatics atom selection
        int noElecStart = options.es1;
        noElecStart = (noElecStart < 1) ? 1 : noElecStart;

        int noElecStop = options.ef1;
        noElecStop = (noElecStop > atoms.length) ? atoms.length : noElecStop;

        for (int i = noElecStart; i <= noElecStop; i++) {
            Atom ai = atoms[i - 1];
            ai.setElectrostatics(false);
            ai.print();
        }

        int nVars = pe.getNumberOfVariables();
        double[] x = new double[nVars];
        pe.getCoordinates(x);

        if (options.gradient) {
            double[] g = new double[nVars];
            int nAts = nVars / 3;
            if (pe instanceof ForceFieldEnergyOpenMM) {
                double[] gOMM = new double[nVars];
                ForceFieldEnergyOpenMM ope = (ForceFieldEnergyOpenMM) pe;
                ope.energyAndGradVsFFX(x, g, gOMM, true);
                for (int i = 0; i < nAts; i++) {
                    int i3 = 3*i;
                    logger.info(String.format(" Atom %d OpenMM gradient:      %14.9g %14.9g %14.9g",
                            i, gOMM[i3], gOMM[i3+1], gOMM[i3+2]));
                    logger.info(String.format(" Atom %d gradient discrepancy: %14.9g %14.9g %14.9g",
                            i, gOMM[i3]-g[i3], gOMM[i3+1]-g[i3+1], gOMM[i3+2]-g[i3+2]));
                }
            } else {
                pe.energyAndGradient(x, g, true);
            }
            for (int i = 0; i < nAts; i++) {
                int i3 = 3*i;
                logger.info(String.format(" Atom %d FFX gradient:         %14.9g %14.9g %14.9g", i, g[i3], g[i3+1], g[i3+2]));
            }
        } else if (pe instanceof ForceFieldEnergyOpenMM) {
            ForceFieldEnergyOpenMM ope = (ForceFieldEnergyOpenMM) pe;
            ope.energyVsFFX(x, true);
        } else {
            pe.energy(x, true);
        }
        //functions.energy(activeAssembly)
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
