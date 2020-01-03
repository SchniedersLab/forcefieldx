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
package ffx.potential.cli;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that save a structure to disc.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class SaveOptions {
    /**
     * -c or --constrain is a flag to print out energy at each step.
     */
    @Option(names = {"-c", "--constrain"}, paramLabel = "false", description = "Apply geometric constraints before saving.")
    private boolean constrain = false;

    private double[] x;
    private double[] outputX;

    /**
     * Performs key operations prior to saving to disc, such as application of geometric constraints.
     *
     * @param mola A MolecularAssembly.
     */
    public void preSaveOperations(MolecularAssembly mola) {
        preSaveOperations(mola.getPotentialEnergy());
    }

    /**
     * Performs key operations prior to saving to disc, such as application of geometric constraints.
     *
     * @param ffe A ForceFieldEnergy.
     */
    public void preSaveOperations(ForceFieldEnergy ffe) {
        if (constrain) {
            int nVars = ffe.getNumberOfVariables();
            if (x == null) {
                x = new double[nVars];
                outputX = new double[nVars];
            }
            x = ffe.getCoordinates(x);
            System.arraycopy(x, 0, outputX, 0, nVars);
            ffe.applyAllConstraintPositions(x, outputX);
            ffe.setCoordinates(outputX);
        }
    }
}
