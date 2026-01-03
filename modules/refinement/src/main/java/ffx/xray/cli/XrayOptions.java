// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
// ******************************************************************************
package ffx.xray.cli;

import ffx.potential.MolecularAssembly;
import ffx.utilities.FFXProperty;
import ffx.xray.solvent.SolventModel;
import ffx.xray.DiffractionData;
import ffx.xray.RefinementEnergy;
import ffx.xray.refine.RefinementMode;
import ffx.xray.parsers.DiffractionFile;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import picocli.CommandLine.ArgGroup;
import picocli.CommandLine.Option;
import picocli.CommandLine.ParseResult;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import static ffx.utilities.PropertyGroup.StructuralRefinement;
import static java.lang.Boolean.parseBoolean;
import static java.lang.Double.parseDouble;
import static java.lang.String.format;

/**
 * Represents command line options for scripts that utilize X-ray data with a maximum likelihood
 * target.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class XrayOptions extends DataRefinementOptions {

  private static final Logger logger = Logger.getLogger(XrayOptions.class.getName());

  /**
   * The ArgGroup keeps the XrayOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n X-ray Refinement Options%n", validate = false)
  private final XrayOptionGroup group = new XrayOptionGroup();
  /**
   * The ArgGroup keeps the XrayReflectionsGroup together when printing help.
   */
  @ArgGroup(heading = "%n X-ray Reflection Data Options%n", validate = false)
  private final XrayReflectionsGroup reflectionGroup = new XrayReflectionsGroup();
  /**
   * The ArgGroup keeps the BFactorGroup together when printing help.
   */
  @ArgGroup(heading = "%n X-ray B-Factor Options%n", validate = false)
  private final BFactorGroup bfactorGroup = new BFactorGroup();
  /**
   * The ArgGroup keeps the BFactorGroup together when printing help.
   */
  @ArgGroup(heading = "%n X-ray Target Options%n", validate = false)
  private final ScatteringGroup targetGroup = new ScatteringGroup();
  /**
   * The ArgGroup keeps the XrayOptionGroup together when printing help.
   */
  @ArgGroup(heading = "%n X-ray Bulk Solvent Options%n", validate = false)
  private final BulkSolventGroup solventGroup = new BulkSolventGroup();
  /**
   * The refinement mode to use.
   */
  public RefinementMode refinementMode = RefinementMode.COORDINATES;
  /**
   * The SolventModel to use.
   */
  public SolventModel solventModel = SolventModel.POLYNOMIAL;

  /**
   * Parse options.
   */
  public void init() {
    refinementMode = RefinementMode.parseMode(group.modeString);
    solventModel = SolventModel.parse(solventGroup.solventString);
  }

  /**
   * Process input to collect Diffraction Files.
   *
   * @param filenames Input filenames (first filename is ignored).
   * @param systems   Currently open systems.
   * @return a list of DiffractionFile instances.
   */
  public List<DiffractionFile> processData(List<String> filenames, MolecularAssembly[] systems) {
    List<DiffractionFile> diffractionfiles = new ArrayList<>();

    logger.info("\n");

    if (filenames.size() > 1) {
      logger.info(format(" Diffraction file = %s, weight = %3.1f, neutron = %b", filenames.get(1), wA, Boolean.FALSE));
      DiffractionFile diffractionfile = new DiffractionFile(filenames.get(1), wA, false);
      diffractionfiles.add(diffractionfile);
    }

    if (reflectionGroup.data != null) {
      for (int i = 0; i < reflectionGroup.data.length; i += 3) {
        boolean neutron = false;
        double w = wA;
        if (reflectionGroup.data.length > i + 1) {
          try {
            w = parseDouble(reflectionGroup.data[i + 1]);
          } catch (Exception e) {
            //
          }
        }
        if (reflectionGroup.data.length > i + 2) {
          try {
            neutron = parseBoolean(reflectionGroup.data[i + 2]);
          } catch (Exception e) {
            //
          }
        }
        logger.info(format(" Diffraction file = %s, weight = %3.1f, neutron = %b", reflectionGroup.data[i], w, neutron));
        DiffractionFile diffractionfile = new DiffractionFile(reflectionGroup.data[i], w, neutron);
        diffractionfiles.add(diffractionfile);
      }

    }

    if (diffractionfiles.isEmpty()) {
      String filename = systems[0].getFile().getAbsolutePath();
      filename = FilenameUtils.removeExtension(filename);
      filename = FilenameUtils.getBaseName(filename);

      logger.info(format(" Diffraction from = %s, weight = %3.1f, neutron = %b", filename, wA, false));
      DiffractionFile diffractionfile = new DiffractionFile(systems, 1.0, false);
      diffractionfiles.add(diffractionfile);
    }

    return diffractionfiles;
  }

  /**
   * setProperties.
   *
   * @param parseResult a {@link picocli.CommandLine.ParseResult} object.
   * @param properties  a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
   */
  public void setProperties(ParseResult parseResult, CompositeConfiguration properties) {
    // wA
    if (!parseResult.hasMatchedOption("wA")) {
      wA = properties.getDouble("data-weight", wA);
    }
    properties.setProperty("data-weight", wA);

    // bSimWeight
    if (!parseResult.hasMatchedOption("bSimWeight")) {
      bfactorGroup.bSimWeight = properties.getDouble("b-sim-weight", bfactorGroup.bSimWeight);
    }
    properties.setProperty("b-sim-weight", bfactorGroup.bSimWeight);

    // F/SigF Cutoff
    if (!parseResult.hasMatchedOption("fSigFCutoff")) {
      reflectionGroup.fSigFCutoff = properties.getDouble("f-sigf-cutoff", reflectionGroup.fSigFCutoff);
    }
    properties.setProperty("f-sigf-cutoff", reflectionGroup.fSigFCutoff);

    // Solvent Grid Search
    if (!parseResult.hasMatchedOption("gridSearch")) {
      solventGroup.gridSearch = properties.getBoolean("solvent-grid-search", solventGroup.gridSearch);
    }
    properties.setProperty("solvent-grid-search", solventGroup.gridSearch);

    // Number of Bins
    if (!parseResult.hasMatchedOption("nBins")) {
      reflectionGroup.nBins = properties.getInt("reflection-bins", reflectionGroup.nBins);
    }
    properties.setProperty("reflection-bins", reflectionGroup.nBins);

    // Grid Sampling
    if (!parseResult.hasMatchedOption("sampling")) {
      targetGroup.sampling = properties.getDouble("sampling", targetGroup.sampling);
    }
    properties.setProperty("sampling", targetGroup.sampling);

    // Atomic Radius Buffer
    if (!parseResult.hasMatchedOption("aRadBuffer")) {
      targetGroup.aRadBuffer = properties.getDouble("scattering-buffer", targetGroup.aRadBuffer);
    }
    properties.setProperty("scattering-buffer", targetGroup.aRadBuffer);

    // RFreeFlag
    if (!parseResult.hasMatchedOption("rFreeFlag")) {
      reflectionGroup.rFreeFlag = properties.getInt("rfree-flag", reflectionGroup.rFreeFlag);
    }
    properties.setProperty("rfree-flag", reflectionGroup.rFreeFlag);

    // Spline Fit
    if (!parseResult.hasMatchedOption("noSplineFit")) {
      targetGroup.noSplineFit = properties.getBoolean("no-spline-fit", targetGroup.noSplineFit);
    }
    properties.setProperty("no-spline-fit", targetGroup.noSplineFit);

    // Use All Gaussians
    if (!parseResult.hasMatchedOption("allGaussians")) {
      targetGroup.useAllGaussians = !properties.getBoolean("use-3g", !targetGroup.useAllGaussians);
    }
    properties.setProperty("use-3g", !targetGroup.useAllGaussians);

    // X-ray Scale Tolerance
    if (!parseResult.hasMatchedOption("xrayScaleTol")) {
      targetGroup.xrayScaleTol = properties.getDouble("xray-scale-tol", targetGroup.xrayScaleTol);
    }
    properties.setProperty("xray-scale-tol", targetGroup.xrayScaleTol);

    // Sigma A Tolerance
    if (!parseResult.hasMatchedOption("sigmaATol")) {
      targetGroup.sigmaATol = properties.getDouble("sigmaa-tol", targetGroup.sigmaATol);
    }
    properties.setProperty("sigmaa-tol", targetGroup.sigmaATol);

    // Number of Residues per B-Factor
    if (!parseResult.hasMatchedOption("nResidueBFactor")) {
      bfactorGroup.nResidueBFactor = properties
          .getInt("n-residue-bfactor", bfactorGroup.nResidueBFactor);
    }
    properties.setProperty("n-residue-bfactor", Integer.toString(bfactorGroup.nResidueBFactor));
    if (bfactorGroup.nResidueBFactor > 0) {
      properties.setProperty("residue-bfactor", "true");
    }

    // Add AnisoU B-Factors to the Refinement
    if (!parseResult.hasMatchedOption("anisoU")) {
      bfactorGroup.anisoU = properties.getBoolean("add-anisou", bfactorGroup.anisoU);
    }
    properties.setProperty("add-anisou", bfactorGroup.anisoU);

    // Refine Molecular Occupancies
    if (!parseResult.hasMatchedOption("refineMolOcc")) {
      group.refineMolOcc = properties.getBoolean("refine-mol-occ", group.refineMolOcc);
    }
    properties.setProperty("refine-mol-occ", group.refineMolOcc);
  }

  /**
   * Process input from opened molecular assemblies to a DiffractionData.
   *
   * @param filenames   All filenames included in the diffraction data.
   * @param assemblies  All molecular assemblies included in the diffraction data.
   * @param properties  The properties to apply.
   * @return An assembled DiffractionData
   */
  public DiffractionData getDiffractionData(List<String> filenames, MolecularAssembly[] assemblies,
      CompositeConfiguration properties) {
    // Set up diffraction data (can be multiple files)
    List<DiffractionFile> diffractionFiles = processData(filenames, assemblies);
    DiffractionData diffractionData = new DiffractionData(assemblies, properties, solventModel,
        diffractionFiles.toArray(new DiffractionFile[0]));
    diffractionData.getRefinementModel().setRefinementMode(refinementMode);
    return diffractionData;
  }

  /**
   * Collection of X-ray Refinement Options.
   */
  private static class XrayOptionGroup {

    /**
     * --rmo or --refineMolOcc
     */
    @Option(names = {"--rmo", "--refineMolOcc"}, paramLabel = "false", defaultValue = "false",
        description = "Refine molecular occupancy.")
    @FFXProperty(name = "refine-mol-occ", propertyGroup = StructuralRefinement, defaultValue = "false",
        description = "Refine molecular occupancy.")
    boolean refineMolOcc = false;

    /**
     * -m or --mode sets the desired refinement mode [COORDINATES, BFACTORS,
     * COORDINATES_AND_BFACTORS, OCCUPANCIES, BFACTORS_AND_OCCUPANCIES, COORDINATES_AND_OCCUPANCIES,
     * COORDINATES_AND_BFACTORS_AND_OCCUPANCIES].
     */
    @Option(names = {"-m", "--mode"}, paramLabel = "coordinates", defaultValue = "coordinates",
        description = "Refinement mode: coordinates, bfactors and/or occupancies.")
    String modeString = "coordinates";
  }

  /**
   * Collection of X-ray Reflection Options.
   */
  private static class XrayReflectionsGroup {

    /**
     * --nBins sets the number of reflection bins to use.
     */
    @Option(names = {"--nBins"}, paramLabel = "10", defaultValue = "10",
        description = "The number of refection bins.")
    @FFXProperty(name = "reflection-bins", propertyGroup = StructuralRefinement, defaultValue = "10",
        description = "The number of refection bins.")
    int nBins = 10;

    /**
     * --FSigFCutoff
     */
    @Option(names = {"--FSigFCutoff"}, paramLabel = "-1.0", defaultValue = "-1.0",
        description = "F / SigF cutoff (-1.0 is no cutoff).")
    @FFXProperty(name = "f-sigf-cutoff", propertyGroup = StructuralRefinement, defaultValue = "-1.0",
        description = "F / SigF cutoff (-1.0 is no cutoff)."
    )
    double fSigFCutoff = -1.0;

    /**
     * -R or --rFreeFlag
     */
    @Option(names = {"-R", "--rFreeFlag"}, paramLabel = "-1", defaultValue = "-1",
        description = "R-Free Flag value (-1 attempts to auto-determine from the data).")
    @FFXProperty(name = "rfree-flag", propertyGroup = StructuralRefinement, defaultValue = "-1",
        description = "R-Free Flag value (-1 attempts to auto-determine from the data)."
    )
    int rFreeFlag = -1;

    /**
     * -X or --data Specify input data filename, weight applied to the data (wA) and if the data is
     * from a neutron experiment.
     */
    @Option(names = {"-X", "--data"}, arity = "3",
        description = "Specify input data filename, its weight (wA) and if its from a neutron experiment (e.g. -X filename 1.0 false).")
    String[] data = null;
  }

  /**
   * Collection of B-Factor Options.
   */
  private static class BFactorGroup {

    /**
     * -B or --bSimWeight
     */
    @Option(names = {"-B", "--bSimWeight"}, paramLabel = "1.0", defaultValue = "1.0",
        description = "B-Factor similarity weight.")
    @FFXProperty(name = "b-sim-weight", propertyGroup = StructuralRefinement, defaultValue = "1.0", description = """
        B-factor harmonic restraint weight between bonded atoms.
        This can be increased for low resolution structures to ~5-10.
        """)
    double bSimWeight = 1.0;

    /**
     * --nResidueBFactor
     */
    @Option(names = {"--nResidueBFactor"}, paramLabel = "0", defaultValue = "0",
        description = "Number of residues per B-factor. 0 uses atomic B-factors (default).")
    @FFXProperty(name = "n-residue-bfactor", propertyGroup = StructuralRefinement, defaultValue = "0",
        description = "Number of residues per B-factor. 0 uses atomic B-factors (default).")
    int nResidueBFactor = 0;

    /**
     * -u or --addAnisoU
     */
    @Option(names = {"-U", "--addAnisoU"}, paramLabel = "false", defaultValue = "false",
        description = "Add Anisotropic B-Factors to refinement.")
    @FFXProperty(name = "add-anisou", propertyGroup = StructuralRefinement, defaultValue = "false",
        description = "Add Anisotropic B-Factors to refinement.")
    boolean anisoU = false;

  }

  /**
   * Collection of X-ray Scattering Options.
   */
  private static class ScatteringGroup {

    /**
     * --aRadBuffer
     */
    @Option(names = {"--aRadBuffer"}, paramLabel = "0.75", defaultValue = "0.75",
        description = "Scattering is evaluated within the atomic radius plus this buffer (Å).")
    @FFXProperty(name = "scattering-buffer", propertyGroup = StructuralRefinement, defaultValue = "0.6",
        description = "Scattering is evaluated within the atomic radius plus this buffer (Å).")
    double aRadBuffer = 0.75;

    /**
     * -G or --sampling
     */
    @Option(names = {"-G", "--sampling"}, paramLabel = "0.6", defaultValue = "0.6",
        description = "The number of grid spaces per Angstrom for the scattering FFT grid.")
    @FFXProperty(name = "sampling", propertyGroup = StructuralRefinement, defaultValue = "0.6",
        description = "The number of grid spaces per Angstrom for the scattering FFT grid.")
    double sampling = 0.6;

    /**
     * --sf or --splineFit
     */
    @Option(names = {"--nsf", "--noSplineFit"}, paramLabel = "false", defaultValue = "false",
        description = "Use a resolution dependent spline scale factor.")
    @FFXProperty(name = "no-spline-fit", propertyGroup = StructuralRefinement, defaultValue = "false",
        description = "Do not use a resolution dependent spline scale factor.")
    boolean noSplineFit = false;

    /**
     * -A or --allGaussians
     */
    @Option(names = {"-A", "--allGaussians"}, paramLabel = "false", defaultValue = "false",
        description = "Use all defined Gaussians for atomic scattering density (the default is to use the top 3).")
    @FFXProperty(name = "use-3g", propertyGroup = StructuralRefinement, defaultValue = "true",
        description = "The three Gaussians with the largest amplitudes define the atomic scattering density.")
    boolean useAllGaussians = false;

    /**
     * --xrayScaleTol
     */
    @Option(names = {"--xrayScaleTol"}, paramLabel = "1.0e-4", defaultValue = "1.0e-4",
        description = "X-ray scale optimization tolerance.")
    @FFXProperty(name = "xray-scale-tol", propertyGroup = StructuralRefinement, defaultValue = "1.0e-4",
        description = "X-ray scale optimization tolerance.")
    double xrayScaleTol = 1.0e-4;

    /**
     * --sigmaATol
     */
    @Option(names = {"--sigmaATol"}, paramLabel = "0.05", defaultValue = "0.05",
        description = "Sigma A optimization tolerance.")
    @FFXProperty(name = "sigmaa-tol", propertyGroup = StructuralRefinement, defaultValue = "0.05",
        description = "Sigma A optimization tolerance.")
    double sigmaATol = 0.05;
  }

  /**
   * Collection of Bulk Solvent Options.
   */
  private static class BulkSolventGroup {

    /**
     * -S or --solventGridSearch
     */
    @Option(names = {"-S", "--solventGridSearch"}, paramLabel = "false", defaultValue = "false",
        description = "Perform a grid search for optimal bulk solvent parameters.")
    @FFXProperty(name = "solvent-grid-search", propertyGroup = StructuralRefinement, defaultValue = "false",
        description = "Perform a grid search for optimal bulk solvent parameters.")
    boolean gridSearch = false;

    /**
     * -sol or --solvent Bulk solvent scattering model [Polynomial/Gaussian/Binary/None].
     */
    @Option(names = {"--sol", "--solvent"}, paramLabel = "POLYNOMIAL", defaultValue = "POLYNOMIAL",
        description = "Bulk solvent scattering model [Polynomial/Gaussian/Binary/None]")
    String solventString = "POLYNOMIAL";
  }

  /**
   * Process input from opened molecular assemblies and diffraction data to a RefinementEnergy.
   *
   * @param diffractionData Diffraction data.
   * @return An assembled RefinementEnergy with X-ray energy.
   */
  public RefinementEnergy toXrayEnergy(DiffractionData diffractionData) {
    diffractionData.scaleBulkFit();
    diffractionData.printStats();
    return new RefinementEnergy(diffractionData);
  }
}
