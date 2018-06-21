package ffx.potential

import ffx.numerics.Potential
import ffx.potential.bonded.Atom
import ffx.potential.bonded.LambdaInterface
import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The CoordShakeEnergy script evaluates the energy of a system before and after moving coordinates by a fixed offset.
 * Intended primarily to test coordinate restraints. Could easily be extended to do randomized offsets, etc.
 * <br>
 * Usage:
 * <br>
 * ffxc CoordShakeEnergy &lt;filename&gt;
 */
@Command(description = " Compute the energy after shaking coordinates.", name = "ffxc CoordShakeEnergy")
class CoordShakeEnergy extends PotentialScript {

    /**
     * -d or --deltaX to set the distance by which the coordinates should be moved. Applies a constant offset in x,
     * reverts, in y, reverts, then z, and reverts again.
     */
    @Option(names = ['-d', '--deltaX'], paramLabel = "1.0", description = 'Distance to move coordinates (applies to all equally)')
    double dX = 1.0
    /**
     * -l or --lambda sets the lambda value to evaluate at.
     */
    @Option(names = ['-l', '--lambda'], paramLabel = "-1", description = 'Lambda value (-1 indicates no lambda dependence).')
    double initialLambda = -1
    /**
     * -s or --start sets the first atom to move.
     */
    @Option(names = ['-s', '--start'], paramLabel = '1', description = 'Starting atom to test')
    int start = 1
    /**
     * -f or --final sets the last atom to move.
     */
    @Option(names = ['-f', '--final'], paramLabel = '-1', description = 'Last atom to test (-1 indicates last atom in the system).')
    int last = -1

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    /**
     * Execute the script.
     */
    def run() {

        if (!init()) {
            return
        }

        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = potentialFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        }

        String modelFilename = activeAssembly.getFile().getAbsolutePath()

        logger.info("\n Running CoordShakeEnergy on " + modelFilename)

        // Turn on computation of lambda derivatives if softcore atoms exist or a single topology has lambda specified.
        boolean lambdaTerm = initialLambda >= 0.0
        if (lambdaTerm) {
            System.setProperty("lambdaterm", "true")
        }
        if (initialLambda < 0.0 || initialLambda > 1.0) {
            initialLambda = 0.0
        }

        Potential thePotential = activeAssembly.getPotentialEnergy()

        if (lambdaTerm) {
            ((LambdaInterface) thePotential).setLambda(initialLambda);
        }

        def eFunct
        if (thePotential instanceof ForceFieldEnergyOpenMM) {
            ForceFieldEnergyOpenMM ommE = (ForceFieldEnergyOpenMM) thePotential
            eFunct = { double[] coords -> return ommE.energyVsFFX(coords, true) }
        } else {
            eFunct = { double[] coords -> return thePotential.energy(coords, true) }
        }

        logger.info(" Starting energy: ")
        double[] x = thePotential.getCoordinates()
        eFunct(x)

        Atom[] atoms = activeAssembly.getAtomArray()

        def axes = ["X", "Y", "Z"]
        double[] atomXYZ = new double[3]
        int firstAt = start - 1
        int lastAt = last
        if (lastAt < 0 || lastAt > atoms.length) {
            lastAt = atoms.length
        }

        for (int i = 0; i < 3; i++) {
            logger.info(String.format(" Moving atoms +1 unit in %s", axes[i]))
            for (int j = firstAt; j < lastAt; j++) {
                Atom atom = atoms[j]
                atom.getXYZ(atomXYZ)
                atomXYZ[i] += dX
                atom.setXYZ(atomXYZ)
            }
            thePotential.getCoordinates(x)
            eFunct(x)

            logger.info(String.format(" Moving atoms -1 unit in %s", axes[i]))
            for (int j = firstAt; j < lastAt; j++) {
                Atom atom = atoms[j]
                atom.getXYZ(atomXYZ)
                atomXYZ[i] -= (2.0 * dX)
                atom.setXYZ(atomXYZ)
            }
            thePotential.getCoordinates(x)
            eFunct(x)

            for (int j = firstAt; j < lastAt; j++) {
                Atom atom = atoms[j]
                atom.getXYZ(atomXYZ)
                atomXYZ[i] += dX
                atom.setXYZ(atomXYZ)
            }
        }

        logger.info(" After movement (may have small rounding errors): ")
        thePotential.getCoordinates(x)
        eFunct(x)
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
