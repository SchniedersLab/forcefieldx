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
package ffx.xray.cli;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import ffx.potential.MolecularAssembly;
import ffx.realspace.RealSpaceFile;
import ffx.xray.RefinementMinimize.RefinementMode;

import picocli.CommandLine.Option;

/**
 * Xray options shared by Xray scripts that use the Pico CLI.
 */
public class RealSpaceOptions {

    private static final Logger logger = Logger.getLogger(RealSpaceOptions.class.getName());

    /**
     * -X or --wA The weight of the real space data (wA).
     */
    @Option(names = {"-X", "--wA"}, paramLabel = "1.0",
            description = "The weight of the real space data (wA).")
    double wA = 1.0;

    /**
     * -x or --data Specify input data filename, weight applied to the data (wA) and if the data is from a neutron experiment.
     */
    @Option(names = {"-x", "--data"}, split = ",",
            description = "Specify input data filename and its weight (wA) (e.g. -X filename,1.0).")
    String[] data = null;

    /**
     * -r or --mode sets the desired refinement mode
     * [COORDINATES, BFACTORS, COORDINATES_AND_BFACTORS, OCCUPANCIES, BFACTORS_AND_OCCUPANCIES, COORDINATES_AND_OCCUPANCIES, COORDINATES_AND_BFACTORS_AND_OCCUPANCIES].
     */
    @Option(names = {"-m", "--mode"}, paramLabel = "coordinates",
            description = "Refinement mode: coordinates, bfactors and/or occupancies.")
    String modeString = "coordinates";

    /**
     * The refinement mode to use.
     */
    RefinementMode refinementMode = RefinementMode.COORDINATES;

    /**
     * Parse options.
     */
    public void init() {
        refinementMode = ffx.xray.RefinementMinimize.parseMode(modeString);
    }

    /**
     * Process input to collect Real Space Files.
     *
     * @param filenames Input filenames (first filename is ignored).
     * @param systems   Currently open systems.
     * @return a list of Real Space File instances.
     */
    public List<RealSpaceFile> processData(List<String> filenames, MolecularAssembly systems[]) {

        logger.info("\n");

        // set up real space map data (can be multiple files)
        List<RealSpaceFile> mapfiles = new ArrayList<>();
        if (filenames.size() > 1) {
            RealSpaceFile realspacefile = new RealSpaceFile(filenames.get(1), wA);
            mapfiles.add(realspacefile);
        }
        if (data != null) {
            for (int i = 0; i < data.length; i += 2) {
                double w = wA;
                if (data.length > i + 1) {
                    try {
                        w = Double.parseDouble(data[i + 1]);
                    } catch (Exception e) {
                        //
                    }
                }
                ffx.realspace.RealSpaceFile realspacefile = new ffx.realspace.RealSpaceFile(data[i], w);
                mapfiles.add(realspacefile);
            }
        }

        if (mapfiles.size() == 0) {
            RealSpaceFile realspacefile = new RealSpaceFile(systems[0], wA);
            mapfiles.add(realspacefile);
        }

        return mapfiles;
    }
}
