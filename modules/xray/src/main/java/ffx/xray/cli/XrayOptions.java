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

import org.apache.commons.io.FilenameUtils;

import ffx.potential.MolecularAssembly;
import ffx.xray.CrystalReciprocalSpace;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.parsers.DiffractionFile;

import picocli.CommandLine.Option;

/**
 * Xray options shared by Xray scripts that use the Pico CLI.
 */
public class XrayOptions {

    private static final Logger logger = Logger.getLogger(XrayOptions.class.getName());

    /**
     * -s or --solvent Bulk solvent scattering model [Polynomial/Gaussian/Binary/None].
     */
    @Option(names = {"-s", "--solvent"}, paramLabel = "solvent", description = "Bulk solvent scattering model [Polynomial/Gaussian/Binary/None]")
    private String solventString = "POLYNOMIAL";

    /**
     * The SolventModel to use.
     */
    SolventModel solventModel = SolventModel.POLYNOMIAL;

    /**
     * -x or --data Specify input data filename, weight applied to the data (wA) and if the data is from a neutron experiment.
     */
    @Option(names = {"-x", "--data"}, split = ",",
            description = "Specify input data filename, its weight (wA) and its from a neutron experiment (e.g. -D filename,1.0,false).")
    String[] data = null;

    /**
     * Parse options.
     */
    public void init() {
        solventModel = CrystalReciprocalSpace.parseSolventModel(solventString);
    }

    /**
     * Process input to collect Diffraction Files.
     *
     * @param filenames Input filenames (first filename is ignored).
     * @param systems   Currently open systems.
     * @return a list of DiffractionFile instances.
     */
    public List<DiffractionFile> processData(List<String> filenames, MolecularAssembly systems[]) {
        List<DiffractionFile> diffractionfiles = new ArrayList<>();

        logger.info("\n");

        if (filenames.size() > 1) {
            logger.info(String.format(" Diffraction file = %s, weight = %3.1f, neutron = %b",
                    filenames.get(1), 1.0, Boolean.FALSE));

            DiffractionFile diffractionfile = new DiffractionFile(filenames.get(1), 1.0, false);
            diffractionfiles.add(diffractionfile);
        }

        if (data != null) {
            for (int i = 0; i < data.length; i += 3) {
                double wA = 1.0;
                boolean neutron = false;
                if (data.length > i + 1) {
                    try {
                        wA = Double.parseDouble(data[i + 1]);
                    } catch (Exception e) {
                        //
                    }
                }
                if (data.length > i + 2) {
                    try {
                        neutron = Boolean.parseBoolean(data[i + 2]);
                    } catch (Exception e) {
                        //
                    }
                }

                logger.info(String.format(" Diffraction file = %s, weight = %3.1f, neutron = %b", data[i], wA, neutron));

                DiffractionFile diffractionfile = new DiffractionFile(data[i], wA, neutron);
                diffractionfiles.add(diffractionfile);
            }
        }

        if (diffractionfiles.size() == 0) {
            String filename = systems[0].getFile().getAbsolutePath();
            filename = FilenameUtils.removeExtension(filename);
            filename = FilenameUtils.getBaseName(filename);

            logger.info(String.format(" Diffraction from = %s, weight = %3.1f, neutron = %b", filename, 1.0, false));
            DiffractionFile diffractionfile = new DiffractionFile(systems, 1.0, false);
            diffractionfiles.add(diffractionfile);
        }

        return diffractionfiles;
    }
}
