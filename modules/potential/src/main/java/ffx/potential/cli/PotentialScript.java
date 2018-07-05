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
package ffx.potential.cli;

import ffx.potential.MolecularAssembly;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;
import ffx.utilities.BaseScript;

/**
 * Base class for scripts in the Potentials package, providing some key functions.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PotentialScript extends BaseScript {

    /**
     * An instance of PotentialFunctions passed into the current context.
     */
    public PotentialsFunctions potentialFunctions;

    /**
     * An active MolecularAssembly passed into the current context or loaded by the Script from a file argument.
     */
    public MolecularAssembly activeAssembly;

    /**
     * Execute the BaseScript init method, then load potential functions.
     *
     * @return Returns true if the script should continue.
     */
    @Override
    public boolean init() {
        if (!super.init()) {
            return false;
        }

        if (context.hasVariable("functions")) {
            potentialFunctions = (PotentialsFunctions) context.getVariable("functions");
        } else {
            potentialFunctions = new PotentialsUtils();
        }

        activeAssembly = null;
        if (context.hasVariable("active")) {
            activeAssembly = (MolecularAssembly) context.getVariable("active");
        }

        return true;
    }

}
