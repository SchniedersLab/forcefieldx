//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.algorithms.optimize.manybody;

import java.util.Arrays;
import static java.lang.System.arraycopy;

import ffx.algorithms.mc.BoltzmannMC;
import ffx.algorithms.optimize.RotamerOptimization;
import ffx.potential.bonded.Residue;

/**
 * Monte Carlo driver for DEE-MC.
 */
public class RotamerMatrixMC extends BoltzmannMC {

    private RotamerOptimization rotamerOptimization;
    private boolean useFullAMOEBAEnergy;
    private final int[] currentRots;
    private final int[] oldRots;
    private final int nRes;
    private final Residue[] residues;


    /**
     * The Rotamers array must be the same array as passed to any MCMove
     * objects used (and not a copy).
     *
     * @param rotamers
     * @param residues
     */
    public RotamerMatrixMC(int[] rotamers, Residue[] residues, boolean useFullAMOEBAEnergy,
                           RotamerOptimization rotamerOptimization) {
        currentRots = rotamers; // This is intentional.
        nRes = rotamers.length;
        oldRots = new int[nRes];
        arraycopy(rotamers, 0, oldRots, 0, nRes);
        this.residues = residues;
        this.useFullAMOEBAEnergy = useFullAMOEBAEnergy;
        this.rotamerOptimization = rotamerOptimization;
    }

    @Override
    public void revertStep() {
        arraycopy(oldRots, 0, currentRots, 0, nRes);
    }

    /**
     * If useFullAMOEBAEnergy is set to true, explicitly evaluates energy,
     * else computes energy from the rotamer energy matrices.
     *
     * @return Energy at the current state
     */
    @Override
    protected double currentEnergy() {
        try {
            try {
                return useFullAMOEBAEnergy ?
                        rotamerOptimization.currentEnergyWrapper(Arrays.asList(residues))
                        : rotamerOptimization.computeEnergy(residues, currentRots, false);
            } catch (ArithmeticException ex) {
                return 1E100;
            }
        } catch (NullPointerException ex) {
            // If using the rotamer energy matrix, and there is some missing
            // energy term, just return a default, very large energy.
            return 1e100;
        }
    }

    @Override
    protected void storeState() {
        arraycopy(currentRots, 0, oldRots, 0, nRes);
    }
}
