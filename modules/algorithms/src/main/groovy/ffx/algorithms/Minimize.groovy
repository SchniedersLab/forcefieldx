package ffx.algorithms

import java.util.stream.Collectors

import org.apache.commons.configuration.CompositeConfiguration
import org.apache.commons.io.FilenameUtils

import edu.rit.pj.ParallelTeam

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.MinimizeOptions
import ffx.numerics.Potential
import ffx.numerics.UnivariateSwitchingFunction
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.LambdaInterface
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.cli.TopologyOptions

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

/**
 * The Minimize script uses a limited-memory BFGS algorithm to minimize the 
 * energy of a molecular system or supersystem.
 * <br>
 * Usage:
 * <br>
 * ffxc Minimizer [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Run L-BFGS minimization on a system.", name = "ffxc Minimize")
class Minimizer extends AlgorithmsScript {

    @Mixin
    MinimizeOptions minimizeOptions

    @Mixin
    AlchemicalOptions alchemical

    @Mixin
    TopologyOptions topology

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = 'Atomic coordinate files in PDB or XYZ format.')
    List<String> filenames = null

    private int threadsAvail = ParallelTeam.getDefaultThreadCount()
    private int threadsPer = threadsAvail
    MolecularAssembly[] topologies


    def run() {

        if (!init()) {
            return
        }

        List<String> arguments = filenames;
        // Check nArgs should either be number of arguments (min 1), else 1.
        int nArgs = arguments ? arguments.size() : 1;
        nArgs = (nArgs < 1) ? 1 : nArgs;

        topologies = new MolecularAssembly[nArgs];

        int numParallel = topology.getNumParallel(threadsAvail, nArgs);
        threadsPer = threadsAvail / numParallel;

        // Turn on computation of lambda derivatives if softcore atoms exist.
        boolean lambdaTerm = alchemical.hasSoftcore() || topology.hasSoftcore();

        if (lambdaTerm) {
            System.setProperty("lambdaterm", "true")
        }

        double lambda = alchemical.getInitialLambda();

        // Relative free energies via the DualTopologyEnergy class require different
        // default OSRW parameters than absolute free energies.
        if (nArgs >= 2) {
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec", "false");
        }

        List<MolecularAssembly> topologyList = new ArrayList<>(4);

        /**
         * Read in files.
         */
        if (!arguments || arguments.isEmpty()) {
            MolecularAssembly mola = algorithmFunctions.getActiveAssembly();
            if (mola == null) {
                return helpString();
            }
            arguments = new ArrayList<>();
            arguments.add(mola.getFile().getName());
            topologyList.add(alchemical.processFile(Optional.of(topology), mola, 0));
        } else {
            logger.info(String.format(" Initializing %d topologies...", nArgs));
            for (int i = 0; i < nArgs; i++) {
                topologyList.add(alchemical.openFile(algorithmFunctions, Optional.of(topology), threadsPer, arguments.get(i), i));
            }
        }

        MolecularAssembly[] topologies = topologyList.toArray(new MolecularAssembly[topologyList.size()]);

        /**
         * Configure the potential to test.
         */
        StringBuilder sb = new StringBuilder("\n Minimizing energy of ");
        Potential potential = topology.assemblePotential(topologies, threadsAvail, sb);

        logger.info(sb.toString());

        LambdaInterface linter = (potential instanceof LambdaInterface) ? (LambdaInterface) potential : null;
        linter?.setLambda(lambda);

        double[] x = new double[potential.getNumberOfVariables()]
        potential.getCoordinates(x)
        potential.energy(x, true)

        AlgorithmListener theListener = algorithmFunctions.getDefaultListener()
        Minimize minimize = new Minimize(topologies[0], potential, theListener)
        minimize.minimize(minimizeOptions.getEps(), minimizeOptions.getIterations());

        potential.getCoordinates(x)
        potential.energy(x, true)

        for (mola in topologies) {
            String filename = mola.getFile().getName()
            String ext = FilenameUtils.getExtension(filename)
            filename = FilenameUtils.removeExtension(filename)

            if (ext.toUpperCase().contains("XYZ")) {
                algorithmFunctions.saveAsXYZ(mola, new File(filename + ".xyz"))
            } else {
                algorithmFunctions.saveAsPDB(mola, new File(filename + ".pdb"))
            }
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
