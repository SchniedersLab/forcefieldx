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
    CompositeConfiguration[] properties
    ForceFieldEnergy[] energies

    private void openFile(String toOpen, int topNum) {
        algorithmFunctions.openAll(toOpen, threadsPer)
        MolecularAssembly mola = algorithmFunctions.getActiveAssembly()
        processFile(mola, topNum)
    }

    private void processFile(MolecularAssembly mola, int topNum) {

        int remainder = (topNum % 2) + 1
        switch (remainder) {
            case 1:
                alchemical.setAlchemicalAtoms(mola)
                alchemical.setUnchargedAtoms(mola)
                break
            case 2:
                topology.setAlchemicalAtoms(mola)
                topology.setUnchargedAtoms(mola)
        }

        // Turn off checks for overlapping atoms, which is expected for lambda=0.
        ForceFieldEnergy energy = mola.getPotentialEnergy()
        energy.getCrystal().setSpecialPositionCutoff(0.0)

        // Save a reference to the topology.
        properties[topNum] = mola.getProperties()
        topologies[topNum] = mola
        energies[topNum] = energy
    }


    def run() {

        if (!init()) {
            return
        }

        List<String> arguments = filenames
        // Check nArgs should either be number of arguments (min 1), else 1.
        int nArgs = arguments ? arguments.size() : 1
        nArgs = (nArgs < 1) ? 1 : nArgs

        topologies = new MolecularAssembly[nArgs]
        properties = new CompositeConfiguration[nArgs]
        energies = new ForceFieldEnergy[nArgs]

        int numParallel = topology.getNumParallel(threadsAvail, nArgs)
        threadsPer = threadsAvail / numParallel

        // Turn on computation of lambda derivatives if softcore atoms exist or a single topology.
        boolean lambdaTerm = alchemical.ligAt1 || topology.ligAt2 || (alchemical.s1 > 0) || (topology.s2 > 0)

        if (lambdaTerm) {
            System.setProperty("lambdaterm", "true")
        }

        // Relative free energies via the DualTopologyEnergy class require different
        // default OSRW parameters than absolute free energies.
        if (nArgs >= 2) {
            // Ligand vapor electrostatics are not calculated. This cancels when the
            // difference between protein and water environments is considered.
            System.setProperty("ligand-vapor-elec", "false")
        }

        /**
         * Read in files.
         */
        if (!arguments || arguments.isEmpty()) {
            MolecularAssembly mola = algorithmFunctions.getActiveAssembly()
            if (mola == null) {
                return helpString()
            }
            arguments = new ArrayList<>()
            arguments.add(mola.getFile().getName())
            processFile(mola)
        } else {
            logger.info(String.format(" Initializing %d topologies...", nArgs))
            for (int i = 0; i < nArgs; i++) {
                openFile(arguments.get(i), i)
            }
        }

        /**
         * Configure the potential to test.
         */
        StringBuilder sb = new StringBuilder("\n Testing lambda derivatives for ")
        Potential potential
        UnivariateSwitchingFunction sf = topology.getSwitchingFunction();
        switch (nArgs) {
            case 1:
                potential = energies[0]
                alchemical.setActiveAtoms(topologies[0]);
                break
            case 2:
            case 4:
            case 8:
                List<Integer> uniqueA
                List<Integer> uniqueB
                if (nArgs >= 4) {
                    uniqueA = topology.getUniqueAtoms(topologies[0], "A");
                    uniqueB = topology.getUniqueAtoms(topologies[1], "B");
                }
                potential = topology.getTopology(topologies, sf, uniqueA, uniqueB, numParallel, sb)
                break
            default:
                logger.severe(" Must have 1, 2, 4, or 8 topologies!")
                break
        }

        ArrayList<MolecularAssembly> topo = new ArrayList<>(Arrays.asList(topologies));
        sb.append(topo.stream().map { t -> t.getFile().getName() }.collect(Collectors.joining(",", "[", "]")))
        logger.info(sb.toString())

        LambdaInterface linter = (LambdaInterface) potential

        // Turn on computation of lambda derivatives if l >= 0 or > 1 argument
        if (lambdaTerm || nArgs > 1) {
            if (alchemical.initialLambda < 0.0 || alchemical.initialLambda > 1.0) {
                alchemical.initialLambda = 0.5
            }
            linter.setLambda(alchemical.initialLambda)
        }

        double[] x = new double[potential.getNumberOfVariables()]
        potential.getCoordinates(x)
        potential.energy(x, true)

        AlgorithmListener theListener = algorithmFunctions.getDefaultListener()
        Minimize minimize = new Minimize(topologies[0], potential, theListener)
        minimize.minimize(minimizeOptions.eps, minimizeOptions.iterations)

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
