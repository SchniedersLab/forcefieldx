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
package ffx.algorithms.mc;

import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Torsion;

/**
 * Represents a random chi[0] spin of the target residue. For use with
 * RosenbluthRotamerMC.
 *
 * @author S. LuCore
 */
public class RosenbluthChi0Move implements MCMove {

    private static final Logger logger = Logger.getLogger(RosenbluthChi0Move.class.getName());

    private final Residue target;
    private final ResidueState origState;
    private final Rotamer newState;
    public final double theta;

    /**
     * <p>Constructor for RosenbluthChi0Move.</p>
     *
     * @param target a {@link ffx.potential.bonded.Residue} object.
     */
    public RosenbluthChi0Move(Residue target) {
        this.target = target;
        AminoAcid3 name = AminoAcid3.valueOf(target.getName());
        origState = target.storeState();
        double chi[] = RotamerLibrary.measureRotamer(target, false);
        theta = ThreadLocalRandom.current().nextDouble(360.0) - 180;
        chi[0] = theta;
        // Need to add sigma values to construct a new Rotamer with these chis.
        double values[] = new double[chi.length * 2];
        for (int i = 0; i < chi.length; i++) {
            int ii = 2 * i;
            values[ii] = chi[i];
            values[ii + 1] = 0.0;
        }
        newState = new Rotamer(name, values);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Performs the move associated with this MCMove. Also updates chi values in
     * associated Torsion objects.
     */
    @Override
    public void move() {
        RotamerLibrary.applyRotamer(target, newState);
        updateTorsions();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Reverts the last applied move() call.
     */
    @Override
    public void revertMove() {
        target.revertState(origState);
        updateTorsions();
    }

    private void updateTorsions() {
        List<ROLS> torsions = target.getTorsionList();
        for (ROLS rols : torsions) {
            ((Torsion) rols).update();
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return String.format("Rosenbluth Rotamer Move:\n   Res:   %s\n   Theta: %3.2f",
                target.toString(), theta);
    }

}
