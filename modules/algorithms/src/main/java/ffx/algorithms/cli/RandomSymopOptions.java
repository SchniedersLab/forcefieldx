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
package ffx.algorithms.cli;

import ffx.potential.MolecularAssembly;

import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that create randomized unit cells.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class RandomSymopOptions {

    /**
     * --rsym or --randomSymOp Apply a random SymOp with translation range -X/2 .. X/2 (0 for random placement in the unit cell, negative for no SymOp).
     */
    @Option(names = {"--rsym", "--randomSymOp"}, paramLabel = "-1.0",
            description = "Apply a random SymOp with translation range -X/2 .. X/2 (0 for random placement in the unit cell, negative for no SymOp)")
    double symScalar = -1.0;

    /**
     * --ruc or --randomUnitCell random unit cell axes will be used achieve the specified density (g/cc) (no default density).
     */
    @Option(names = {"--ruc", "--randomUnitCell"}, paramLabel = "-1.0",
            description = "Apply random unit cell axes to achieve the specified density (g/cc).")
    double ucDensity = -1.0;

    /**
     * <p>randomize.</p>
     *
     * @param assembly a {@link ffx.potential.MolecularAssembly} object.
     */
    public void randomize(MolecularAssembly assembly) {
        if (symScalar >= 0) {
            assembly.applyRandomSymOp(symScalar);
        }
        if (ucDensity > 0) {
            assembly.applyRandomDensity(ucDensity);
        }
    }
}
