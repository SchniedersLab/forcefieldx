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

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import ffx.potential.MolecularAssembly;
import ffx.xray.CrystalReciprocalSpace;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.parsers.DiffractionFile;
import ffx.xray.RefinementMinimize;
import ffx.xray.RefinementMinimize.RefinementMode;

import picocli.CommandLine.Option;

/**
 * Xray options shared by Xray scripts that use the Pico CLI.
 */
public class XrayOptions {

    private static final Logger logger = Logger.getLogger(XrayOptions.class.getName());

    /**
     * -s or --solvent Bulk solvent scattering model [Polynomial/Gaussian/Binary/None].
     */

    @Option(names = {"-s", "--solvent"}, paramLabel = "POLYNOMIAL",
            description = "Bulk solvent scattering model [Polynomial/Gaussian/Binary/None]")
    private String solventString = "POLYNOMIAL";

    /**
     * The SolventModel to use.
     */
    SolventModel solventModel = SolventModel.POLYNOMIAL;

    /**
     * -S or --solventGridSearch
     */
    @Option(names = {"-S", "--solventGridSearch"}, paramLabel = "false",
            description = "Perform a grid search for optimal bulk solvent parameters.")
    boolean gridSearch = false;

    /**
     * -N or --nBins
     */
    @Option(names = {"-N", "--nBins"}, paramLabel = "10",
            description = "The number of refection bins.")
    int nBins = 10;

    /**
     * -F or --FSigFCutoff
     */
    @Option(names = {"-F", "--FSigFCutoff"}, paramLabel = "-1.0",
            description = "F / SigF cutoff (-1.0 is no cutoff).")
    double fSigFCutoff = -1.0;

    /**
     * -X or --wA The weight of the X-ray data (wA).
     */
    @Option(names = {"-X", "--wA"}, paramLabel = "1.0",
            description = "The weight of the X-ray data (wA).")
    double wA = 1.0;

    /**
     * -a or -aRadBuffer
     */
    @Option(names = {"-a", "--aRadBuffer"}, paramLabel = "0.75",
            description = "Set the distance beyond the atomic radius to evaluate scattering (A).")
    double aRadBuffer = 0.75;

    /**
     * -G or --sampling
     */
    @Option(names = {"-G", "--sampling"}, paramLabel = "0.6", description = "The number of grid spaces per Angstrom for the scattering FFT grid.")
    double sampling = 0.6;

    /**
     * -R or --rFreeFlag
     */
    @Option(names = {"-R", "--rFreeFlag"}, paramLabel = "-1",
            description = "R-Free Flag value (-1 attempts to auto-determine from the data).")
    int rFreeFlag = -1;

    /**
     * --sf or --splineFit
     */
    @Option(names = {"--sf", "--splineFit"}, paramLabel = "true",
            description = "Use a resolution dependent spline scale.")
    boolean splineFit = true;

    /**
     * -A or --all
     */
    @Option(names = {"-A", "--all"}, paramLabel = "false",
            description = "Use all defined Gaussians for atomic scattering density (the default is to use the top 3).")
    boolean useAll = false;

    /**
     * --xrayScaleTol
     */
    @Option(names = {"--xrayScaleTol"}, paramLabel = "1.0e-4", description = "X-ray scale optimization tolerance.")
    double xrayScaleTol = 1.0e-4;

    /**
     * --sigmaATol
     */
    @Option(names = {"--sigmaATol"}, paramLabel = "0.05", description = "Sigma A optimization tolerance.")
    double sigmaATol = 1.0e-4;

    /**
     * -B or --bSimWeight
     */
    @Option(names = {"-B", "--bSimWeight"}, paramLabel = "1.0", description = "B-Factor similarity weight.")
    double bSimWeight = 1.0;

    /**
     * --residueBFactor
     */
    @Option(names = {"--residueBFactor"}, paramLabel = "false", description = "Refine B-factors by residue.")
    boolean residueBFactor = false;

    /**
     * --nResidueBFactor
     */
    @Option(names = {"--nResidueBFactor"}, paramLabel = "1",
            description = "Number of residues per B-factor, if refining by B-factors by residue.")
    int nResidueBFactor = 1;

    /**
     * -u or --addAnisoU
     */
    @Option(names = {"-U", "--addAnisoU"}, paramLabel = "false", description = "Add Anisotropic B-Factors to refinement.")
    boolean anisoU = false;

    /**
     * -o or --refineMolOcc
     */
    @Option(names = {"-o", "--refineMolOcc"}, paramLabel = "false", description = "Refine on molecules.")
    boolean refineMolOcc = false;

    /**
     * -x or --data Specify input data filename, weight applied to the data (wA) and if the data is from a neutron experiment.
     */
    @Option(names = {"-x", "--data"}, split = ",",
            description = "Specify input data filename, its weight (wA) and if its from a neutron experiment (e.g. -D filename,1.0,false).")
    String[] data = null;

    /**
     * -r or --mode sets the desired refinement mode
     * [COORDINATES, BFACTORS, COORDINATES_AND_BFACTORS, OCCUPANCIES, BFACTORS_AND_OCCUPANCIES, COORDINATES_AND_OCCUPANCIES, COORDINATES_AND_BFACTORS_AND_OCCUPANCIES].
     */
    @Option(names = {"-m", "--mode"}, paramLabel = "coordinates", description = "Refinement mode: coordinates, bfactors and/or occupancies.")
    String modeString = "coordinates";

    /**
     * The refinement mode to use.
     */
    RefinementMode refinementMode = RefinementMode.COORDINATES;

    /**
     * Parse options.
     */
    public void init() {
        solventModel = CrystalReciprocalSpace.parseSolventModel(solventString);
        refinementMode = RefinementMinimize.parseMode(modeString);
    }

    public void setProperties(CompositeConfiguration properties) {
        properties.setProperty("xweight", wA);
        properties.setProperty("bsimweight", bSimWeight);
        properties.setProperty("fsigfcutoff", fSigFCutoff);
        properties.setProperty("gridsearch", gridSearch);
        properties.setProperty("nbins", nBins);
        properties.setProperty("sampling", sampling);
        properties.setProperty("aradbuff", aRadBuffer);
        properties.setProperty("rrfreeflag", rFreeFlag);
        properties.setProperty("splinefit", splineFit);
        properties.setProperty("use_3g", !useAll);
        properties.setProperty("xrayscaletol", xrayScaleTol);
        properties.setProperty("sigmaatol", sigmaATol);
        properties.setProperty("residuebfactor", residueBFactor);
        properties.setProperty("nresiduebfactor", nResidueBFactor);
        properties.setProperty("addanisou", anisoU);
        properties.setProperty("refinemolocc", refineMolOcc);
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
                    filenames.get(1), wA, Boolean.FALSE));

            DiffractionFile diffractionfile = new DiffractionFile(filenames.get(1), wA, false);
            diffractionfiles.add(diffractionfile);
        }

        if (data != null) {
            for (int i = 0; i < data.length; i += 3) {
                boolean neutron = false;
                double w = wA;
                if (data.length > i + 1) {
                    try {
                        w = Double.parseDouble(data[i + 1]);
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
                logger.info(String.format(" Diffraction file = %s, weight = %3.1f, neutron = %b", data[i], w, neutron));
                DiffractionFile diffractionfile = new DiffractionFile(data[i], w, neutron);
                diffractionfiles.add(diffractionfile);
            }
        }

        if (diffractionfiles.size() == 0) {
            String filename = systems[0].getFile().getAbsolutePath();
            filename = FilenameUtils.removeExtension(filename);
            filename = FilenameUtils.getBaseName(filename);

            logger.info(String.format(" Diffraction from = %s, weight = %3.1f, neutron = %b", filename, wA, false));
            DiffractionFile diffractionfile = new DiffractionFile(systems, 1.0, false);
            diffractionfiles.add(diffractionfile);
        }

        return diffractionfiles;
    }
}
