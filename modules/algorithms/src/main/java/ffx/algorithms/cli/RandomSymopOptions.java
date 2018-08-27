/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.algorithms.cli;

import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.crystal.SymOp;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;

import picocli.CommandLine;

/**
 * Represents command line options for scripts that create randomized unit cells.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class RandomSymopOptions {

    private static final Logger logger = Logger.getLogger(RandomSymopOptions.class.getName());

    /**
     * --rsym or --randomSymOp to apply a random Cartesian symmetry operator with the specified translation range -X .. X (no default).
     */
    @CommandLine.Option(names = {"--rsym", "--randomSymOp"}, paramLabel = "-1.0",
            description = "Apply a random Cartesian symmetry operator with a random translation in the range -X .. X; < 0 disables.")
    double symScalar = -1.0;

    /**
     * --ruc or --randomUnitCell random unit cell axes will be used achieve the specified density (g/cc) (no default density).
     */
    @CommandLine.Option(names = {"--ruc", "--randomUnitCell"}, paramLabel = "-1.0",
            description = "Apply random unit cell axes to achieve the specified density (g/cc).")
    double ucDensity = -1.0;

    /**
     * <p>randomize.</p>
     *
     * @param assembly  a {@link ffx.potential.MolecularAssembly} object.
     * @param potential a {@link ffx.crystal.CrystalPotential} object.
     */
    public void randomize(MolecularAssembly assembly, CrystalPotential potential) {
        applyRandomSymop(assembly);
        applyRandomDensity(assembly, potential);
    }

    /**
     * Randomizes position in the unit cell by applying a Cartesian symop with a random translation.
     *
     * @param assembly Assembly to randomize position of.
     */
    public void applyRandomSymop(MolecularAssembly assembly) {
        if (symScalar > 0) {
            SymOp symOp = SymOp.randomSymOpFactory(symScalar);
            logger.info(String.format("\n Applying random Cartesian SymOp:\n%s", symOp.toString()));
            Crystal crystal = assembly.getCrystal();
            Atom[] atoms = assembly.getAtomArray();
            double[] xyz = new double[3];
            for (int i = 0; i < atoms.length; i++) {
                atoms[i].getXYZ(xyz);
                crystal.applyCartesianSymOp(xyz, xyz, symOp);
                atoms[i].setXYZ(xyz);
            }
        }
    }

    /**
     * Applies a randomly drawn density to a molecular system's crystal.
     *
     * @param assembly Assembly to randomize the density of.
     * @throws java.lang.IllegalArgumentException if any.
     */
    public void applyRandomDensity(MolecularAssembly assembly) throws IllegalArgumentException {
        applyRandomDensity(assembly, assembly.getPotentialEnergy());
    }

    /**
     * Applies a randomly drawn density to a molecular system's crystal.
     *
     * @param assembly  Assembly to randomize the density of.
     * @param potential CrystalPotential to apply the new Crystal parameters to.
     */
    public void applyRandomDensity(MolecularAssembly assembly, CrystalPotential potential) {
        if (ucDensity > 0) {
            logger.info(String.format("\n Applying random unit cell axes with target density of %6.3f\n",
                    ucDensity));
            Crystal crystal = assembly.getCrystal();
            if (!crystal.aperiodic()) {
                double mass = assembly.getMass();
                crystal.randomParameters(ucDensity, mass);
                potential.setCrystal(crystal);
            } else {
                logger.fine(String.format(" Potential %s is an aperiodic system!", potential));
            }
        }
    }
}
