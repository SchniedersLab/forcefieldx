// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.xray;

import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import ffx.xray.parallel.GridMethod;
import ffx.xray.scatter.NeutronFormFactor;
import ffx.xray.solvent.SolventModel;
import ffx.xray.parsers.CCP4MapWriter;
import ffx.xray.parsers.DiffractionFile;
import ffx.xray.parsers.MTZWriter;
import ffx.xray.parsers.MTZWriter.MTZType;
import ffx.xray.refine.RefinementMode;
import ffx.xray.refine.RefinementModel;
import ffx.xray.scatter.XRayFormFactor;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.numerics.math.ScalarMath.b2u;
import static ffx.utilities.TinkerUtils.version;
import static java.lang.Math.abs;
import static java.lang.String.format;
import static java.util.Arrays.fill;
import static org.apache.commons.io.FilenameUtils.removeExtension;

/**
 * DiffractionData class.
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class DiffractionData implements DataContainer {

  private static final Logger logger = Logger.getLogger(DiffractionData.class.getName());
  private final MolecularAssembly[] assembly;
  private final String modelName;
  private final DiffractionFile[] dataFiles;
  private final int n;
  private final Crystal[] crystal;
  private final Resolution[] resolution;
  private final ReflectionList[] reflectionList;
  private final DiffractionRefinementData[] refinementData;
  private final CrystalReciprocalSpace[] crystalReciprocalSpacesFc;
  private final CrystalReciprocalSpace[] crystalReciprocalSpacesFs;
  private final SolventModel solventModel;
  private final RefinementModel refinementModel;
  private final boolean use_3g;
  private final double aRadBuff;
  private final double xrayScaleTol;
  private final double sigmaATol;
  private final double bSimWeight;
  private final double bNonZeroWeight;
  private final boolean nativeEnvironmentApproximation;
  private final SigmaAMinimize[] sigmaAMinimize;
  private final CrystalStats[] crystalStats;
  private ParallelTeam parallelTeam;
  private final GridMethod gridMethod;
  private final boolean[] scaled;
  private double xWeight;
  /**
   * If true, perform a grid search for bulk solvent parameters.
   */
  private final boolean gridSearch;
  /**
   * If true, fit a scaling spline between Fo and Fc.
   */
  private final boolean splineFit;

  /**
   * construct a diffraction data assembly
   *
   * @param molecularAssemblies {@link ffx.potential.MolecularAssembly molecular assembly} object array
   *                            (typically containing alternate conformer assemblies), used as the atomic model for
   *                            comparison against the data
   * @param properties          system properties file
   * @param solventModel        the type of solvent model desired - see {@link
   *                            SolventModel bulk solvent model} selections
   * @param dataFiles            one or more {@link ffx.xray.parsers.DiffractionFile} to be refined against
   */
  public DiffractionData(MolecularAssembly[] molecularAssemblies, CompositeConfiguration properties,
      SolventModel solventModel, DiffractionFile... dataFiles) {
    this.assembly = molecularAssemblies;
    this.solventModel = solventModel;
    this.modelName = molecularAssemblies[0].getName();
    this.dataFiles = dataFiles;
    this.n = dataFiles.length;

    int rflag = properties.getInt("rfree-flag", -1);
    // Settings
    double fsigfCutoff = properties.getDouble("f-sigf-cutoff", -1.0);
    gridSearch = properties.getBoolean("solvent-grid-search", false);
    splineFit = !properties.getBoolean("no-spline-fit", false);
    use_3g = properties.getBoolean("use-3g", true);
    aRadBuff = properties.getDouble("scattering-buffer", 0.75);
    double sampling = properties.getDouble("sampling", 0.6);
    xrayScaleTol = properties.getDouble("xray-scale-tol", 1e-4);
    sigmaATol = properties.getDouble("sigmaa-tol", 0.05);
    xWeight = properties.getDouble("data-weight", 1.0);
    bSimWeight = properties.getDouble("b-sim-weight", 1.0);
    bNonZeroWeight = properties.getDouble("b-nonzero-weight", 1.0);
    boolean residueBFactor = properties.getBoolean("residue-bfactor", false);
    int nResidueBFactor = properties.getInt("n-residue-bfactor", 1);
    boolean addAnisou = properties.getBoolean("add-anisou", false);
    boolean refineMolOcc = properties.getBoolean("refine-mol-occ", false);

    ForceField forceField = molecularAssemblies[0].getForceField();
    nativeEnvironmentApproximation =
        forceField.getBoolean("NATIVE_ENVIRONMENT_APPROXIMATION", false);

    crystal = new Crystal[n];
    resolution = new Resolution[n];
    reflectionList = new ReflectionList[n];
    refinementData = new DiffractionRefinementData[n];
    sigmaAMinimize = new SigmaAMinimize[n];
    crystalStats = new CrystalStats[n];

    for (int i = 0; i < n; i++) {
      crystal[i] = Crystal.checkProperties(properties);
      if (crystal[i] == null) {
        crystal[i] = molecularAssemblies[i].getCrystal().getUnitCell();
      }
      File dataFile = new File(dataFiles[i].getFilename());
      double res = dataFiles[i].getDiffractionfilter().getResolution(dataFile, crystal[i]);
      if (dataFiles[i].isNeutron()) {
        resolution[i] = Resolution.checkProperties(properties, true, res);
      } else {
        resolution[i] = Resolution.checkProperties(properties, false, res);
      }
      if (resolution[i] == null) {
        logger.severe(" Resolution could not be determined from property or reflection files.");
      }
      reflectionList[i] = new ReflectionList(crystal[i], resolution[i], properties);

      logger.info(resolution[i].toString());

      refinementData[i] = new DiffractionRefinementData(properties, reflectionList[i]);
      dataFiles[i].getDiffractionfilter().readFile(dataFile, reflectionList[i], refinementData[i], properties);
    }

    if (logger.isLoggable(Level.INFO)) {
      StringBuilder sb = new StringBuilder();
      sb.append("\n X-ray Refinement Settings\n\n");
      sb.append("  Target Function\n");
      sb.append("   X-ray refinement weight: ").append(xWeight).append("\n");
      sb.append("   Use cctbx 3 Gaussians: ").append(use_3g).append("\n");
      sb.append("   Atomic form factor radius buffer: ").append(aRadBuff).append("\n");
      sb.append("   Reciprocal space sampling rate: ").append(sampling).append("\n");
      sb.append("   Resolution dependent spline scale: ").append(splineFit).append("\n");
      sb.append("   Solvent grid search: ").append(gridSearch).append("\n");
      sb.append("   X-ray scale fit tolerance: ").append(xrayScaleTol).append("\n");
      sb.append("   Sigma A fit tolerance: ").append(sigmaATol).append("\n");
      sb.append("   Native environment approximation: ").append(nativeEnvironmentApproximation)
          .append("\n");
      sb.append("  Reflections\n");
      sb.append("   F/sigF cutoff: ").append(fsigfCutoff).append("\n");
      sb.append("   R Free flag (-1 auto-determine from the data): ").append(rflag).append("\n");
      sb.append("   Number of bins: ").append(reflectionList[0].nBins).append("\n");
      sb.append("  B-Factors\n");
      sb.append("   Similarity weight: ").append(bSimWeight).append("\n");
      // sb.append("  Non-zero weight (bnonzeroweight): ").append(bNonZeroWeight).append("\n");
      // sb.append("  Lagrangian mass (bmass): ").append(bMass).append("\n");
      sb.append("   Refine by residue: ").append(residueBFactor).append("\n");
      if (residueBFactor) {
        sb.append("   Number of residues per B: ").append(nResidueBFactor).append(")\n");
      }
      sb.append("   Add ANISOU for refinement: ").append(addAnisou).append("\n");
      sb.append("  Occupancies\n");
      sb.append("   Refine on molecules: ").append(refineMolOcc).append("\n");
      // sb.append("  Lagrangian mass (occmass): ").append(occMass).append("\n");

      logger.info(sb.toString());
    }

    // now set up the refinement model
    refinementModel = new RefinementModel(molecularAssemblies);

    // initialize atomic form factors
    // int index = 0;
    for (Atom a : refinementModel.getScatteringAtoms()) {
      // logger.info(" Diffraction Atom " + index + ": " + a);
      // index++;

      XRayFormFactor atomff = new XRayFormFactor(a, use_3g, 2.0);
      // NeutronFormFactor neutronFF = new NeutronFormFactor(a, 2.0);
      if (a.getOccupancy() == 0.0) {
        a.setFormFactorWidth(1.0);
        continue;
      }

      double arad;
      try {
        arad = a.getVDWType().radius * 0.5;
      } catch (NullPointerException ex) {
        logger.warning(format(" Failure to get van der Waals type for atom %s; ensure the vdW term is enabled!", a));
        throw ex;
      }
      double[] xyz = new double[3];
      xyz[0] = a.getX() + arad;
      xyz[1] = a.getY();
      xyz[2] = a.getZ();
      double rho = atomff.rho(0.0, 1.0, xyz);
      // double rhoN = abs(neutronFF.rho(0.0, 1.0, xyz));
      // while (rho > 0.001 || rhoN > 0.001) {
      while (rho > 0.001) {
        arad += 0.1;
        xyz[0] = a.getX() + arad;
        rho = atomff.rho(0.0, 1.0, xyz);
        // rhoN = abs(neutronFF.rho(0.0, 1.0, xyz));
      }
      arad += aRadBuff;
      a.setFormFactorWidth(arad);
    }

    // set up FFT and run it
    crystalReciprocalSpacesFc = new CrystalReciprocalSpace[n];
    crystalReciprocalSpacesFs = new CrystalReciprocalSpace[n];

    parallelTeam = molecularAssemblies[0].getParallelTeam();

    String gridString = properties.getString("grid-method", "SLICE");
    gridMethod = GridMethod.parse(gridString);

    parallelTeam = new ParallelTeam();
    for (int i = 0; i < n; i++) {
      // Atomic Scattering
      crystalReciprocalSpacesFc[i] = new CrystalReciprocalSpace(
              reflectionList[i], refinementModel.getScatteringAtoms(),
              parallelTeam, parallelTeam, false,
              this.dataFiles[i].isNeutron(), SolventModel.NONE, gridMethod);
      refinementData[i].setCrystalReciprocalSpaceFc(crystalReciprocalSpacesFc[i]);
      crystalReciprocalSpacesFc[i].setUse3G(use_3g);
      crystalReciprocalSpacesFc[i].setWeight(this.dataFiles[i].getWeight());
      crystalReciprocalSpacesFc[i].lambdaTerm = false;
      crystalReciprocalSpacesFc[i].setNativeEnvironmentApproximation(
          nativeEnvironmentApproximation);

      // Bulk Solvent Scattering
      crystalReciprocalSpacesFs[i] = new CrystalReciprocalSpace(
              reflectionList[i], refinementModel.getScatteringAtoms(),
              parallelTeam, parallelTeam, true,
              this.dataFiles[i].isNeutron(), solventModel, gridMethod);
      refinementData[i].setCrystalReciprocalSpaceFs(crystalReciprocalSpacesFs[i]);
      crystalReciprocalSpacesFs[i].setUse3G(use_3g);
      crystalReciprocalSpacesFs[i].setWeight(this.dataFiles[i].getWeight());
      crystalReciprocalSpacesFs[i].lambdaTerm = false;
      crystalReciprocalSpacesFs[i].setNativeEnvironmentApproximation(
          nativeEnvironmentApproximation);
      crystalStats[i] = new CrystalStats(reflectionList[i], refinementData[i]);
    }

    scaled = new boolean[n];
    fill(scaled, false);
  }

  /**
   * read in a different assembly to average in structure factors
   *
   * @param assembly the {@link ffx.potential.MolecularAssembly molecular assembly} object array
   *                 (typically containing alternate conformer assemblies), used as the atomic model to average
   *                 in with previous assembly
   * @param index    the current data index (for cumulative average purposes)
   */
  public void AverageFc(MolecularAssembly[] assembly, int index) {
    RefinementModel tmprefinementmodel = new RefinementModel(assembly);

    // initialize atomic form factors
    for (Atom a : tmprefinementmodel.getScatteringAtoms()) {
      XRayFormFactor atomff = new XRayFormFactor(a, use_3g, 2.0);
      if (a.getOccupancy() == 0.0) {
        a.setFormFactorWidth(1.0);
        continue;
      }

      double arad = a.getVDWType().radius * 0.5;
      double[] xyz = new double[3];
      xyz[0] = a.getX() + arad;
      xyz[1] = a.getY();
      xyz[2] = a.getZ();
      double rho = atomff.rho(0.0, 1.0, xyz);
      while (rho > 0.001) {
        arad += 0.1;
        xyz[0] = a.getX() + arad;
        rho = atomff.rho(0.0, 1.0, xyz);
      }
      arad += aRadBuff;
      a.setFormFactorWidth(arad);
    }

    // set up FFT and run it
    for (int i = 0; i < n; i++) {
      crystalReciprocalSpacesFc[i] =
          new CrystalReciprocalSpace(
              reflectionList[i],
              tmprefinementmodel.getScatteringAtoms(),
              parallelTeam,
              parallelTeam,
              false,
              dataFiles[i].isNeutron(),
              SolventModel.NONE,
              gridMethod);
      crystalReciprocalSpacesFc[i].setNativeEnvironmentApproximation(
          nativeEnvironmentApproximation);
      refinementData[i].setCrystalReciprocalSpaceFc(crystalReciprocalSpacesFc[i]);

      crystalReciprocalSpacesFs[i] =
          new CrystalReciprocalSpace(
              reflectionList[i],
              tmprefinementmodel.getScatteringAtoms(),
              parallelTeam,
              parallelTeam,
              true,
              dataFiles[i].isNeutron(),
              solventModel,
              gridMethod);
      crystalReciprocalSpacesFs[i].setNativeEnvironmentApproximation(
          nativeEnvironmentApproximation);
      refinementData[i].setCrystalReciprocalSpaceFs(crystalReciprocalSpacesFs[i]);
    }

    int nhkl = refinementData[0].n;
    double[][] fc = new double[nhkl][2];
    double[][] fs = new double[nhkl][2];

    // Run FFTs.
    for (int i = 0; i < n; i++) {
      crystalReciprocalSpacesFc[i].computeDensity(fc);
      if (solventModel != SolventModel.NONE) {
        crystalReciprocalSpacesFs[i].computeDensity(fs);
      }

      // Average in with current data.
      for (int j = 0; j < refinementData[i].n; j++) {
        refinementData[i].fc[j][0] += (fc[j][0] - refinementData[i].fc[j][0]) / index;
        refinementData[i].fc[j][1] += (fc[j][1] - refinementData[i].fc[j][1]) / index;

        refinementData[i].fs[j][0] += (fs[j][0] - refinementData[i].fs[j][0]) / index;
        refinementData[i].fs[j][1] += (fs[j][1] - refinementData[i].fs[j][1]) / index;
      }
    }
  }

  /**
   * Parallelized call to compute atomic density on a grid, followed by FFT to compute structure
   * factors.
   *
   * @see CrystalReciprocalSpace#computeDensity(double[][], boolean)
   */
  public void computeAtomicDensity() {
    for (int i = 0; i < n; i++) {
      crystalReciprocalSpacesFc[i].computeDensity(refinementData[i].fc);
      if (solventModel != SolventModel.NONE) {
        crystalReciprocalSpacesFs[i].computeDensity(refinementData[i].fs);
      }
    }
  }

  /**
   * Similar to Potential.destroy(), frees up resources associated with this RealSpaceData.
   *
   * @return If assets successfully freed.
   */
  public boolean destroy() {
    try {
      boolean underlyingShutdown = true;
      for (MolecularAssembly assem : assembly) {
        // Continue trying to shut assemblies down even if one fails to shut down.
        boolean thisShutdown = assem.destroy();
        underlyingShutdown = underlyingShutdown && thisShutdown;
      }
      parallelTeam.shutdown();
      return underlyingShutdown;
    } catch (Exception ex) {
      logger.warning(format(" Exception in shutting down a RealSpaceData: %s", ex));
      logger.info(Utilities.stackTraceToString(ex));
      return false;
    }
  }

  /**
   * Getter for the field <code>assembly</code>.
   *
   * @return the assembly
   */
  public MolecularAssembly[] getAssembly() {
    return assembly;
  }

  /**
   * Getter for the field <code>crystal</code>.
   *
   * @return the crystal
   */
  public Crystal[] getCrystal() {
    return crystal;
  }

  /**
   * Getter for the field <code>crs_fc</code>.
   *
   * @return the crs_fc
   */
  public CrystalReciprocalSpace[] getCrystalReciprocalSpacesFc() {
    return crystalReciprocalSpacesFc;
  }

  /**
   * Getter for the field <code>crs_fs</code>.
   *
   * @return the crs_fs
   */
  public CrystalReciprocalSpace[] getCrystalReciprocalSpacesFs() {
    return crystalReciprocalSpacesFs;
  }

  /**
   * Getter for the field <code>n</code>.
   *
   * @return the n
   */
  public int getN() {
    return n;
  }

  /**
   * Getter for the field <code>parallelTeam</code>.
   *
   * @return the parallelTeam
   */
  public ParallelTeam getParallelTeam() {
    return parallelTeam;
  }

  /**
   * Getter for the field <code>refinementData</code>.
   *
   * @return the refinementData
   */
  public DiffractionRefinementData[] getRefinementData() {
    return refinementData;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public RefinementModel getRefinementModel() {
    return refinementModel;
  }

  /**
   * Getter for the field <code>reflectionList</code>.
   *
   * @return the reflectionList
   */
  public ReflectionList[] getReflectionList() {
    return reflectionList;
  }

  /**
   * Getter for the field <code>resolution</code>.
   *
   * @return the resolution
   */
  public Resolution[] getResolution() {
    return resolution;
  }

  /**
   * Getter for the field <code>scaled</code>.
   *
   * @return the scaled
   */
  public boolean[] getScaled() {
    return scaled;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getWeight() {
    return xWeight;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setWeight(double weight) {
    this.xWeight = weight;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String printEnergyUpdate() {
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < n; i++) {
      sb.append(format(
              "     dataset %d (weight: %5.1f): R: %6.2f Rfree: %6.2f chemical energy: %8.2f likelihood: %8.2f free likelihood: %8.2f\n",
              i + 1,
              dataFiles[i].getWeight(),
              crystalStats[i].getR(),
              crystalStats[i].getRFree(),
              assembly[0].getPotentialEnergy().getTotalEnergy(),
              dataFiles[i].getWeight() * sigmaAMinimize[i].calculateLikelihood(),
              dataFiles[i].getWeight() * refinementData[i].llkF));
    }

    return sb.toString();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String printOptimizationHeader() {
    return "R  Rfree";
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String printOptimizationUpdate() {
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < n; i++) {
      sb.append(format("%6.2f %6.2f ", crystalStats[i].getR(), crystalStats[i].getRFree()));
    }

    return sb.toString();
  }

  /**
   * Print all statistics for all datasets associated with the model.
   */
  public void printStats() {
    int numberOfScatteringAtoms = 0;
    int nHeavyScatteringAtoms = 0;

    // Note - we are including inactive and/or un-used atoms in the following loop.
    for (Atom a : refinementModel.getScatteringAtoms()) {
      if (a.getOccupancy() == 0.0) {
        continue;
      }
      numberOfScatteringAtoms++;
      if (a.getAtomicNumber() == 1) {
        continue;
      }
      nHeavyScatteringAtoms++;
    }

    for (int i = 0; i < n; i++) {
      if (!scaled[i]) {
        scaleBulkFit(i);
      }

      String sb = format(" Statistics for Data Set %d of %d\n\n"
                  + "  Weight:     %6.2f\n  Neutron data: %4s\n"
                  + "  Model:        %s\n  Data file:    %s\n",
              i + 1, n, dataFiles[i].getWeight(), dataFiles[i].isNeutron(),
              modelName, dataFiles[i].getFilename());
      logger.info(sb);

      crystalStats[i].printScaleStats();
      crystalStats[i].printDPIStats(nHeavyScatteringAtoms, numberOfScatteringAtoms);
      crystalStats[i].printHKLStats();
      crystalStats[i].printSNStats();
      crystalStats[i].printRStats();
    }
  }

  /**
   * Scale model and fit bulk solvent to all data.
   */
  public void scaleBulkFit() {
    for (int i = 0; i < n; i++) {
      scaleBulkFit(i);
    }
  }

  /**
   * Scale model and fit bulk solvent to dataset i of n.
   *
   * @param i a int.
   */
  public void scaleBulkFit(int i) {
    String sb = format(
            " Scaling Data Set %d of %d\n\n  Weight: %6.2f\n  Neutron data: %s\n  Model: %s\n  Data file: %s",
            i + 1, n, dataFiles[i].getWeight(), dataFiles[i].isNeutron(),
            modelName, dataFiles[i].getFilename());
    logger.info(sb);


    // Reset Overall Scale Factor K.
    refinementData[i].modelScaleK = 0.0;
    // Reset Overall B-Factor.
    Arrays.fill(refinementData[i].modelAnisoB, 0.0);

    // Reset Bulk Solvent K and B.
    if (dataFiles[i].isNeutron()) {
      // Water Deuterium scattering
      // Deuterium Water Neutron Scattering 5.803 + 2 * 6.671 = 19.145
      // Hydrogen Water Neutron Scattering 5.803 - 2 * -3.7390 = -1.675
      refinementData[i].bulkSolventK = 0.639;
    } else {
      // 8 oxygen electrons and 2 hydrogen electrons per water = 10
      // Electron density of water in units of electrons / A^3
      refinementData[i].bulkSolventK = 0.334;
    }
    refinementData[i].bulkSolventUeq = b2u(50.0);

    // Set Bulk Solvent Model A and B Parameters.
    refinementData[i].crystalReciprocalSpaceFs.setDefaultSolventAB();

    // Check for Fixed Bulk Solvent K and B.
    CompositeConfiguration properties = assembly[0].getProperties();
    if (dataFiles[i].isNeutron()) {
      if (properties.containsKey("neutron-bulk-solvent")) {
        String string = properties.getString("neutron-bulk-solvent", "0.639 50.0");
        String[] split = string.split("\\s+");
        refinementData[i].bulkSolventK = Double.parseDouble(split[0]);
        refinementData[i].bulkSolventUeq = b2u(Double.parseDouble(split[1]));
        refinementData[i].bulkSolventFixed = true;
      } else {
        refinementData[i].bulkSolventFixed = false;
      }
    } else {
      if (properties.containsKey("bulk-solvent")) {
        String string = properties.getString("bulk-solvent", "0.334 50.0");
        String[] split = string.split("\\s+");
        refinementData[i].bulkSolventK = Double.parseDouble(split[0]);
        refinementData[i].bulkSolventUeq = b2u(Double.parseDouble(split[1]));
        refinementData[i].bulkSolventFixed = true;
      } else {
        refinementData[i].bulkSolventFixed = false;
      }
    }

    // Update Computed Structure Factors.
    crystalReciprocalSpacesFc[i].computeDensity(refinementData[i].fc);
    if (solventModel != SolventModel.NONE) {
      crystalReciprocalSpacesFs[i].computeDensity(refinementData[i].fs);
    }

    // Optimize Bulk Solvent Model.
    if (solventModel != SolventModel.NONE) {
      boolean bulkSolventFixed = refinementData[i].bulkSolventFixed;

      // First, optimize the overall model K and B factor with bulk solvent fixed.
      refinementData[i].bulkSolventFixed = true;
      ScaleBulkMinimize scaleBulkMinimize = new ScaleBulkMinimize(reflectionList[i],
          refinementData[i], crystalReciprocalSpacesFs[i], parallelTeam);
      scaleBulkMinimize.minimize(xrayScaleTol);

      // Bulk solvent grid searches.
      if (!bulkSolventFixed) {
        refinementData[i].bulkSolventFixed = false;
        scaleBulkMinimize = new ScaleBulkMinimize(reflectionList[i],
            refinementData[i], crystalReciprocalSpacesFs[i], parallelTeam);
      }

      // The search for bulk solvent model parameters (e.g., probe radius and shrink radius)
      // can be slow due to recomputing the bulk solvent scattering factors.
      if (gridSearch) {
        // Search for bulk solvent model parameters.
        scaleBulkMinimize.gridOptimizeBulkSolventModel();
      }

      // Search for overall Ks and Bs.
      if (!bulkSolventFixed) {
        scaleBulkMinimize.gridOptimizeKsBs();
      }

      // Final optimization of Overall K/B and Bulk Solvent K/B.
      scaleBulkMinimize.minimize(xrayScaleTol);
    } else {
      // Optimization of both overall and bulk solvent parameters.
      ScaleBulkMinimize scaleBulkMinimize = new ScaleBulkMinimize(
          reflectionList[i], refinementData[i], crystalReciprocalSpacesFs[i], parallelTeam);
      scaleBulkMinimize.minimize(xrayScaleTol);
    }

    // sigmaA / LLK calculation
    sigmaAMinimize[i] = new SigmaAMinimize(reflectionList[i], refinementData[i], parallelTeam);
    sigmaAMinimize[i].minimize(sigmaATol);

    if (splineFit) {
      // initialize minimizers
      SplineMinimize splineMinimize = new SplineMinimize(
          reflectionList[i], refinementData[i], refinementData[i].spline, SplineEnergy.SplineType.FOFC);
      splineMinimize.minimize(1e-5);
    }

    scaled[i] = true;
  }

  /**
   * Perform 10 Fc calculations for the purposes of timings.
   */
  public void timings() {
    logger.info(" Performing 10 Fc calculations for timing...");
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < 10; j++) {
        crystalReciprocalSpacesFc[i].computeDensity(refinementData[i].fc, true);
        crystalReciprocalSpacesFs[i].computeDensity(refinementData[i].fs, true);
        crystalReciprocalSpacesFc[i].computeAtomicGradients(
            refinementData[i].dFc, refinementData[i].freeR, refinementData[i].rFreeFlag,
            RefinementMode.COORDINATES, true);
        crystalReciprocalSpacesFs[i].computeAtomicGradients(
            refinementData[i].dFs, refinementData[i].freeR, refinementData[i].rFreeFlag,
            RefinementMode.COORDINATES, true);
      }
    }
  }

  /**
   * write current datasets to MTZ files
   *
   * @param filename output filename, or filename root for multiple datasets
   * @see MTZWriter#write()
   */
  public void writeData(String filename) {
    if (n == 1) {
      writeData(filename, 0);
    } else {
      for (int i = 0; i < n; i++) {
        writeData(removeExtension(filename) + "_" + i + ".mtz", i);
      }
    }
  }

  /**
   * write dataset i to MTZ file
   *
   * @param filename output filename
   * @param i        dataset to write out
   * @see MTZWriter#write()
   */
  public void writeData(String filename, int i) {
    MTZWriter mtzwriter;
    if (scaled[i]) {
      mtzwriter = new MTZWriter(reflectionList[i], refinementData[i], filename);
    } else {
      mtzwriter = new MTZWriter(reflectionList[i], refinementData[i], filename, MTZType.DATAONLY);
    }
    mtzwriter.write();
  }

  /**
   * write 2Fo-Fc and Fo-Fc maps for all datasets
   *
   * @param filename output root filename for Fo-Fc and 2Fo-Fc maps
   */
  public void writeMaps(String filename, boolean normalize) {
    if (n == 1) {
      writeMaps(filename, 0, normalize);
    } else {
      for (int i = 0; i < n; i++) {
        writeMaps(removeExtension(filename) + "_" + i + ".map", i, normalize);
      }
    }
  }

  /**
   * write 2Fo-Fc and Fo-Fc maps for a datasets
   *
   * @param filename output root filename for Fo-Fc and 2Fo-Fc maps
   * @param i        a int.
   * @param normalize Normalize the maps.
   */
  private void writeMaps(String filename, int i, boolean normalize) {
    if (!scaled[i]) {
      scaleBulkFit(i);
    }

    // Fo-Fc
    crystalReciprocalSpacesFc[i].computeAtomicGradients(refinementData[i].foFc1,
        refinementData[i].freeR, refinementData[i].rFreeFlag, RefinementMode.COORDINATES);
    double[] densityGrid = crystalReciprocalSpacesFc[i].getDensityGrid();
    int extx = (int) crystalReciprocalSpacesFc[i].getXDim();
    int exty = (int) crystalReciprocalSpacesFc[i].getYDim();
    int extz = (int) crystalReciprocalSpacesFc[i].getZDim();

    CCP4MapWriter mapwriter = new CCP4MapWriter(extx, exty, extz, crystal[i],
        removeExtension(filename) + "_fofc.map");
    mapwriter.write(densityGrid, normalize);

    // 2Fo-Fc
    crystalReciprocalSpacesFc[i].computeAtomicGradients(refinementData[i].foFc2,
        refinementData[i].freeR, refinementData[i].rFreeFlag, RefinementMode.COORDINATES);
    densityGrid = crystalReciprocalSpacesFc[i].getDensityGrid();
    extx = (int) crystalReciprocalSpacesFc[i].getXDim();
    exty = (int) crystalReciprocalSpacesFc[i].getYDim();
    extz = (int) crystalReciprocalSpacesFc[i].getZDim();

    mapwriter = new CCP4MapWriter(extx, exty, extz, crystal[i],
        removeExtension(filename) + "_2fofc.map");
    mapwriter.write(densityGrid, normalize);
  }

  /**
   * Write current model to PDB file.
   *
   * @param filename output PDB filename
   */
  public void writeModel(String filename) {
    StringBuilder remark = new StringBuilder();

    File file = version(new File(filename));
    CompositeConfiguration properties = assembly[0].getProperties();
    ForceField forceField = assembly[0].getForceField();
    PDBFilter pdbFilter = new PDBFilter(file, Arrays.asList(assembly), forceField, properties);

    Date now = new Date();
    SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss ");
    remark.append("REMARK FFX output ISO-8601 date: ").append(sdf.format(now)).append("\n");
    remark.append("REMARK\n");
    remark.append("REMARK   3\n");
    remark.append("REMARK   3 REFINEMENT\n");
    remark.append("REMARK   3   PROGRAM     : FORCE FIELD X\n");
    remark.append("REMARK   3\n");
    for (int i = 0; i < n; i++) {
      remark.append("REMARK   3  DATA SET ").append(i + 1).append("\n");
      if (dataFiles[i].isNeutron()) {
        remark.append("REMARK   3   DATA SET TYPE   : NEUTRON\n");
      } else {
        remark.append("REMARK   3   DATA SET TYPE   : X-RAY\n");
      }
      remark.append("REMARK   3   DATA SET WEIGHT : ").append(dataFiles[i].getWeight()).append("\n");
      remark.append("REMARK   3\n");
      remark.append(crystalStats[i].getPDBHeaderString());
    }
    for (int i = 0; i < assembly.length; i++) {
      remark.append("REMARK   3  CHEMICAL SYSTEM ").append(i + 1).append("\n");
      remark.append(assembly[i].getPotentialEnergy().getPDBHeaderString());
    }
    pdbFilter.writeFileWithHeader(file, remark);
  }

  /**
   * Write current model to PDB file.
   *
   * @param filename output PDB filename
   */
  public void writeModel(String filename, Set<Atom> excludeAtoms, double pH) {
    StringBuilder remark = new StringBuilder();

    File file = version(new File(filename));
    PDBFilter pdbFilter = new PDBFilter(file, Arrays.asList(assembly), null, null);
    if (pH > 0) {
      pdbFilter.setRotamerTitration(true);
    }

    Date now = new Date();
    SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss ");
    remark.append("REMARK FFX output ISO-8601 date: ").append(sdf.format(now)).append("\n");
    remark.append("REMARK\n");
    remark.append("REMARK   3\n");
    remark.append("REMARK   3 REFINEMENT\n");
    remark.append("REMARK   3   PROGRAM     : FORCE FIELD X\n");
    remark.append("REMARK   3\n");
    for (int i = 0; i < n; i++) {
      remark.append("REMARK   3  DATA SET ").append(i + 1).append("\n");
      if (dataFiles[i].isNeutron()) {
        remark.append("REMARK   3   DATA SET TYPE   : NEUTRON\n");
      } else {
        remark.append("REMARK   3   DATA SET TYPE   : X-RAY\n");
      }
      remark.append("REMARK   3   DATA SET WEIGHT : ").append(dataFiles[i].getWeight()).append("\n");
      remark.append("REMARK   3\n");
      remark.append(crystalStats[i].getPDBHeaderString());
    }
    for (int i = 0; i < assembly.length; i++) {
      remark.append("REMARK   3  CHEMICAL SYSTEM ").append(i + 1).append("\n");
      remark.append(assembly[i].getPotentialEnergy().getPDBHeaderString());
    }
    remark.append("REMARK   3   TITRATION PH   : \n").append(pH).append("\n");
    String[] remarks = remark.toString().split("\n");
    pdbFilter.writeFile(file, false, excludeAtoms, true, true, remarks);
  }

  /**
   * Write bulk solvent mask for all datasets to a CCP4 map file.
   *
   * @param filename output filename, or output root filename for multiple datasets
   * @see CCP4MapWriter#write(double[], boolean)
   */
  public void writeSolventMask(String filename) {
    if (n == 1) {
      writeSolventMask(filename, 0);
    } else {
      for (int i = 0; i < n; i++) {
        writeSolventMask(removeExtension(filename) + "_" + i + ".map", i);
      }
    }
  }

  /**
   * Write bulk solvent mask for all datasets to a CNS map file.
   *
   * @param filename output filename, or output root filename for multiple datasets
   */
  public void writeSolventMaskCNS(String filename) {
    if (n == 1) {
      writeSolventMaskCNS(filename, 0);
    } else {
      for (int i = 0; i < n; i++) {
        writeSolventMaskCNS(removeExtension(filename) + "_" + i + ".map", i);
      }
    }
  }

  /**
   * Move the atomic coordinates for the FFT calculation - this is independent of the {@link
   * ffx.potential.bonded.Atom} object coordinates as it uses a linearized array
   *
   * @see CrystalReciprocalSpace#updateCoordinates()
   */
  void updateCoordinates() {
    for (int i = 0; i < n; i++) {
      crystalReciprocalSpacesFc[i].updateCoordinates();
      crystalReciprocalSpacesFs[i].updateCoordinates();
    }
  }

  /**
   * Determine the total atomic gradients due to the diffraction data - performs an inverse FFT
   *
   * @param refinementMode the {@link RefinementMode refinement mode}
   *                       requested
   * @see CrystalReciprocalSpace#computeAtomicGradients(double[][], int[], int,
   * RefinementMode, boolean) computeAtomicGradients
   */
  void computeAtomicGradients(RefinementMode refinementMode) {
    for (int i = 0; i < n; i++) {
      crystalReciprocalSpacesFc[i].computeAtomicGradients(refinementData[i].dFc,
          refinementData[i].freeR, refinementData[i].rFreeFlag, refinementMode);
      crystalReciprocalSpacesFs[i].computeAtomicGradients(refinementData[i].dFs,
          refinementData[i].freeR, refinementData[i].rFreeFlag, refinementMode);
    }
  }

  /**
   * compute total crystallographic likelihood target
   *
   * @return the total -log likelihood
   * @see SigmaAMinimize#calculateLikelihood()
   */
  double computeLikelihood() {
    double e = 0.0;
    for (int i = 0; i < n; i++) {
      e += dataFiles[i].getWeight() * sigmaAMinimize[i].calculateLikelihood();
    }
    return e;
  }

  /**
   * setLambdaTerm.
   *
   * @param lambdaTerm a boolean.
   */
  void setLambdaTerm(boolean lambdaTerm) {
    for (int i = 0; i < n; i++) {
      crystalReciprocalSpacesFc[i].setLambdaTerm(lambdaTerm);
      crystalReciprocalSpacesFs[i].setLambdaTerm(lambdaTerm);
    }
  }

  /**
   * Write bulk solvent mask for dataset i to a CNS map file.
   *
   * @param filename output filename
   * @param i        dataset to write out
   */
  private void writeSolventMaskCNS(String filename, int i) {
    if (solventModel != SolventModel.NONE) {
      try {
        PrintWriter cnsfile = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
        cnsfile.println(" ANOMalous=FALSE");
        cnsfile.println(" DECLare NAME=FS DOMAin=RECIprocal TYPE=COMP END");
        for (HKL ih : reflectionList[i].hklList) {
          int j = ih.getIndex();
          cnsfile.printf(" INDE %d %d %d FS= %.4f %.4f\n",
              ih.getH(), ih.getK(), ih.getL(),
              refinementData[i].fsF(j), Math.toDegrees(refinementData[i].fsPhi(j)));
        }
        cnsfile.close();
      } catch (Exception e) {
        String message = "Fatal exception evaluating structure factors.\n";
        logger.log(Level.SEVERE, message, e);
        System.exit(-1);
      }
    }
  }

  /**
   * Write bulk solvent mask for dataset i to a CCP4 map file.
   *
   * @param filename output filename
   * @param i        dataset to write out
   * @see CCP4MapWriter#write(double[], boolean)
   */
  private void writeSolventMask(String filename, int i) {
    if (solventModel != SolventModel.NONE) {
      CCP4MapWriter mapwriter =
          new CCP4MapWriter(
              (int) crystalReciprocalSpacesFs[i].getXDim(),
              (int) crystalReciprocalSpacesFs[i].getYDim(),
              (int) crystalReciprocalSpacesFs[i].getZDim(),
              crystal[i], filename);
      mapwriter.write(crystalReciprocalSpacesFs[i].getSolventGrid());
    }
  }

  /**
   * Getter for the field <code>bSimWeight</code>.
   *
   * @return the bSimWeight
   */
  double getbSimWeight() {
    return bSimWeight;
  }

  /**
   * Getter for the field <code>bNonZeroWeight</code>.
   *
   * @return the bNonZeroWeight
   */
  double getbNonZeroWeight() {
    return bNonZeroWeight;
  }

}
