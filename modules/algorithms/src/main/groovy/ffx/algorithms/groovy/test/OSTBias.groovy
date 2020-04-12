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

import org.apache.commons.io.FilenameUtils

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder

import ffx.algorithms.cli.OSTOptions
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering
import ffx.numerics.switching.PowerSwitch
import ffx.potential.DualTopologyEnergy
import ffx.potential.ForceFieldEnergy
import ffx.potential.bonded.Atom

/**
 * The OSTBias script tests the Transition-Tempered Orthogonal Space Random Walk Potential.
 * <br>
 * Usage:
 * <br>
 * ffxc test.OSTBias [options] &lt;filename [file2...]&gt;
 */
class OSTBias extends Script {

    /**
     * Options for the OSTBias Script.
     * <br>
     * Usage:
     * <br>
     * ffxc test.OSTBias [options] &lt;filename [file2...]&gt;
     */
    class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.')
        boolean help;
        /**
         * -s1 or --start1 defines the first softcored atom for the first topology.
         */
        @Option(shortName = 's1', longName = 'start1', defaultValue = '0', description = 'Starting ligand atom for 1st topology')
        int s1;
        /**
         * -s2 or --start2 defines the first softcored atom for the second topology.
         */
        @Option(shortName = 's2', longName = 'start2', defaultValue = '0', description = 'Starting ligand atom for 2nd topology')
        int s2;
        /**
         * -f1 or --final1 defines the last softcored atom for the first topology.
         */
        @Option(shortName = 'f1', longName = 'final1', defaultValue = '-1', description = 'Final ligand atom for the 1st topology')
        int f1;
        /**
         * -f2 or --final2 defines the last softcored atom for the second topology.
         */
        @Option(shortName = 'f2', longName = 'final2', defaultValue = '-1', description = 'Final ligand atom for the 2nd topology')
        int f2;
        /**
         * -es1 or --noElecStart1 defines the first atom of the first topology to have no electrostatics.
         */
        @Option(shortName = 'es1', longName = 'noElecStart1', defaultValue = '1', description = 'Starting no-electrostatics atom for 1st topology')
        int es1;
        /**
         * -es2 or --noElecStart2 defines the first atom of the second topology to have no electrostatics.
         */
        @Option(shortName = 'es2', longName = 'noElecStart2', defaultValue = '1', description = 'Starting no-electrostatics atom for 2nd topology')
        int es2;
        /**
         * -ef1 or --noElecFinal1 defines the last atom of the first topology to have no electrostatics.
         */
        @Option(shortName = 'ef1', longName = 'noElecFinal1', defaultValue = '-1', description = 'Final no-electrostatics atom for 1st topology')
        int ef1;
        /**
         * -ef2 or --noElecFinal2 defines the last atom of the second topology to have no electrostatics.
         */
        @Option(shortName = 'ef2', longName = 'noElecFinal2', defaultValue = '-1', description = 'Final no-electrostatics atom for 2nd topology')
        int ef2;
        /**
         * -l or --lambda sets the lambda value to minimize at.
         */
        @Option(shortName = 'l', longName = 'lambda', defaultValue = '0.5',
                description = 'Initial lambda value.')
        double lambda;
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed
        List<String> filenames;
    }

    @Override
    OSTBias run() {
        def cli = new CliBuilder(usage: ' ffxc test.OSTBias [options] <filename> [file2...]', header: ' Options:');

        def options = new Options();
        cli.parseFromInstance(options, args);

        if (options.help == true) {
            cli.usage()
            return this
        }

        List<String> arguments = options.filenames;

        // Read in command line file.
        String filename = arguments.get(0);

        // First atom of the ligand.
        int ligandStart = options.s1

        // Last atom of the ligand.
        int ligandStop = options.f1

        // First atom of the ligand for the 2nd topology.
        int ligandStart2 = options.s2

        // Last atom of the ligand for the 2nd topology.
        int ligandStop2 = options.f2

        // First atom for no electrostatics.
        int noElecStart = options.es1

        // Last atom for no electrostatics.
        int noElecStop = options.ef1

        // First atom of the 2nd topology for no electrostatics.
        int noElecStart2 = options.es2

        // Last atom of the 2nd topology for no electrostatics.
        int noElecStop2 = options.ef2

        // Initial lambda value (0 is ligand in vacuum; 1 is ligand in PBC).
        double initialLambda = options.lambda;

        println("\n Testing Orthogonal Space Random Walk on " + filename);

        File structureFile = new File(FilenameUtils.normalize(filename));
        structureFile = new File(structureFile.getAbsolutePath());
        String baseFilename = FilenameUtils.removeExtension(structureFile.getName());
        File histogramRestart = new File(baseFilename + ".his");
        File lambdaRestart = null;
        File dyn = null;

        if (!histogramRestart.exists()) {
            logger.severe("\n Histogram restart file does not exist.");
        } else if (!histogramRestart.canRead()) {
            logger.severe("\n Histogram restart file can not be read.");
        }

        // For a single process job, try to get the restart files from the current directory.
        lambdaRestart = new File(baseFilename + ".lam");
        dyn = new File(baseFilename + ".dyn");

        if (!dyn.exists()) {
            dyn = null;
        }

        // Turn on computation of lambda derivatives.
        System.setProperty("lambdaterm", "true");

        // Relative free energies via the DualTopologyEnergy class require different
        // default OST parameters than absolute free energies.
        if (arguments.size() == 2) {
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec", "false");
        }

        // Open the first system
        open(filename);

        // Get a reference to the first system's ForceFieldEnergy and atom array.
        ForceFieldEnergy energy = active.getPotentialEnergy();
        Atom[] atoms = active.getAtomArray();
        // Apply the ligand atom selection
        for (int i = ligandStart; i <= ligandStop; i++) {
            Atom ai = atoms[i - 1];
            ai.setApplyLambda(true);
            ai.print();
        }

        // Apply the no electrostatics atom selection
        if (noElecStart < 1) {
            noElecStart = 1;
        }
        if (noElecStop > atoms.length) {
            noElecStop = atoms.length;
        }
        for (int i = noElecStart; i <= noElecStop; i++) {
            Atom ai = atoms[i - 1];
            ai.setElectrostatics(false);
            ai.print();
        }

        // Turn off checks for overlapping atoms, which is expected for lambda=0.
        energy.getCrystal().setSpecialPositionCutoff(0.0);
        // OST will be configured for either single or dual topology.
        OrthogonalSpaceTempering orthogonalSpaceTempering = null;
        // Save a reference to the first topology.
        topology1 = active;

        if (arguments.size() == 1) {
            orthogonalSpaceTempering = OSTOptions.constructOST(energy, lambdaRestart, histogramRestart, active, null, sh);
        } else {
            // Open the 2nd topology.
            filename = arguments.get(1);
            open(filename);
            energy = active.getPotentialEnergy();
            atoms = active.getAtomArray();
            // Apply the ligand atom selection for the 2nd topology.
            for (int i = ligandStart2; i <= ligandStop2; i++) {
                Atom ai = atoms[i - 1];
                ai.setApplyLambda(true);
                ai.print();
            }

            // Apply the no electrostatics atom selection
            if (noElecStart2 < 1) {
                noElecStart2 = 1;
            }
            if (noElecStop2 > atoms.length) {
                noElecStop2 = atoms.length;
            }
            for (int i = noElecStart2; i <= noElecStop2; i++) {
                Atom ai = atoms[i - 1];
                ai.setElectrostatics(false);
                ai.print();
            }

            // Save a reference to the second topology.
            topology2 = active;
            // Turn off checks for overlapping atoms, which is expected for lambda=0.
            energy.getCrystal().setSpecialPositionCutoff(0.0);
            // Create the DualTopology potential energy.
            DualTopologyEnergy dualTopologyEnergy = new DualTopologyEnergy(topology1, active, new PowerSwitch());
            orthogonalSpaceTempering = OSTOptions.constructOST(energy, lambdaRestart, histogramRestart, active, null, sh);
        }

        /**
         * Stop propagating lambda to prevent adding new Gaussians
         * to the biasing potential, which would introduce artifacts into the
         * finite-difference derivatives.
         */
        orthogonalSpaceTempering.setPropagateLambda(false);
        orthogonalSpaceTempering.setLambda(initialLambda);
        n = orthogonalSpaceTempering.getNumberOfVariables();

        assert (n % 3 == 0);
        n = n / 3;

        // Finite-difference step size.
        double step = 1.0e-5;

        double[] x = new double[3 * n];
        double[] analytic = new double[3 * n];
        double[] g = new double[3 * n];
        double[] numeric = new double[3];
        orthogonalSpaceTempering.getCoordinates(x);

        // Test Lambda gradients.
        for (int j = 0; j < 3; j++) {
            double lambda = initialLambda - 0.001 + 0.001 * j;

            if (lambda - step < 0.0) {
                continue;
            } else if (lambda + step > 1.0) {
                continue;
            } else {
                orthogonalSpaceTempering.setLambda(lambda);
            }

            // Calculate the energy and analytic dE/dX
            double eL = orthogonalSpaceTempering.energyAndGradient(x, g);

            // Analytic dEdL
            double dEdLambda = orthogonalSpaceTempering.getTotaldEdLambda();

            // Calculate the finite-difference dEdL
            orthogonalSpaceTempering.setLambda(lambda + step);
            double lp = orthogonalSpaceTempering.energyAndGradient(x, g);

            orthogonalSpaceTempering.setLambda(lambda - step);
            double lm = orthogonalSpaceTempering.energyAndGradient(x, g);

            double dedl = (lp - lm) / (2.0 * step);

            logger.info(String.format(" Analytic dE/dL:   %15.8f", dEdLambda));
            logger.info(String.format(" Numeric  dE/dL:   %15.8f\n", dedl));

            // Calculate analytic dE/dX/dL
            orthogonalSpaceTempering.setLambda(lambda);
            double e = 0.0;
            double orig = 0.0;
            double gradientTolerance = 1.0e-3;

            // Calculate finite-difference coordinate gradient
            for (int i = 0; i < n; i++) {
                //Atom a0 = atoms[i];
                int i3 = i * 3;
                int i0 = i3 + 0;
                int i1 = i3 + 1;
                int i2 = i3 + 2;

                // Calculate the analytic dE/dX
                orthogonalSpaceTempering.energyAndGradient(x, analytic);

                // Find numeric dX
                orig = x[i0];
                x[i0] = orig + step;
                e = orthogonalSpaceTempering.energyAndGradient(x, g);
                x[i0] = orig - step;
                e = e - orthogonalSpaceTempering.energyAndGradient(x, g);
                x[i0] = orig;
                numeric[0] = e / (2.0 * step);

                // Find numeric dY
                orig = x[i1];
                x[i1] = orig + step;
                e = orthogonalSpaceTempering.energyAndGradient(x, g);
                x[i1] = orig - step;
                e = e - orthogonalSpaceTempering.energyAndGradient(x, g);
                x[i1] = orig;
                numeric[1] = e / (2.0 * step);

                // Find numeric dZ
                orig = x[i2];
                x[i2] = orig + step;
                e = orthogonalSpaceTempering.energyAndGradient(x, g);
                x[i2] = orig - step;
                e = e - orthogonalSpaceTempering.energyAndGradient(x, g);
                x[i2] = orig;
                numeric[2] = e / (2.0 * step);

                double dx = analytic[i0] - numeric[0];
                double dy = analytic[i1] - numeric[1];
                double dz = analytic[i2] - numeric[2];
                double len = Math.sqrt(dx * dx + dy * dy + dz * dz);
                if (len > gradientTolerance) {
                    logger.info(" Atom " + i + String.format(" failed: %10.6f.", len)
                            + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[i0], analytic[i1], analytic[i2])
                            + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]));
                    return;
                } else {
                    logger.info(" Atom " + i + String.format(" passed: %10.6f.", len)
                            + String.format("\n Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[i0], analytic[i1], analytic[i2])
                            + String.format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]));
                }
            }
        }

        return this
    }
}
