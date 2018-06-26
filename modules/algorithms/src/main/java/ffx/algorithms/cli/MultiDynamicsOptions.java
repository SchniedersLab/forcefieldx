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
package ffx.algorithms.cli;

import picocli.CommandLine;

/**
 * Represents command line options for scripts that can create multiple walkers,
 * such as multi-walker OSRW. Should be kept agnostic to whether it is an MD-based
 * algorithm, or some other flavor of Monte Carlo.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class MultiDynamicsOptions {

    /**
     * -S or --synchronous sets synchronous walker communication (not recommended)
     */
    @CommandLine.Option(names = {"-S", "--synchronous"},
            description = "Walker communication is synchronous")
    private boolean synchronous = false;

    /**
     * -dw or --distributeWalkers allows walkers to start from multiple
     * conformations; AUTO picks up per-walker conformations as
     * filename.pdb_(walker number), and specifying a residue starts a
     * rotamer optimization to generate side-chain configurations to start
     * from.
     */
    @CommandLine.Option(names = {"--dw", "--distributeWalkers"}, paramLabel = "OFF", description = "AUTO: Pick up per-walker configurations as [filename.pdb]_[num], or specify a residue to distribute on.")
    private String distributeWalkersString;

    public boolean isSynchronous() {
        return synchronous;
    }
}
