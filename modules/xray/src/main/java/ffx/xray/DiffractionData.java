// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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

import static ffx.xray.CrystalReciprocalSpace.SolventModel.POLYNOMIAL;
import static java.lang.String.format;
import static java.util.Arrays.fill;

import edu.rit.pj.ParallelTeam;
import ffx.crystal.CCP4MapWriter;
import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.RefinementMinimize.RefinementMode;
import ffx.xray.parsers.DiffractionFile;
import ffx.xray.parsers.MTZWriter;
import ffx.xray.parsers.MTZWriter.MTZType;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

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
  // Settings
  private final double fsigfCutoff;
  private final boolean use_3g;
  private final double aRadBuff;
  private final double xrayScaleTol;
  private final double sigmaATol;
  private final double bSimWeight;
  private final double bNonZeroWeight;
  private final double bMass;
  private final boolean residueBFactor;
  private final int nResidueBFactor;
  private final boolean addAnisou;
  private final boolean refineMolOcc;
  private final double occMass;
  private final boolean nativeEnvironmentApproximation;
  private ScaleBulkMinimize[] scaleBulkMinimize;
  private SigmaAMinimize[] sigmaAMinimize;
  private SplineMinimize[] splineMinimize;
  private CrystalStats[] crystalStats;
  private ParallelTeam parallelTeam;
  private CrystalReciprocalSpace.GridMethod gridMethod;
  private boolean[] scaled;
  private double xWeight;
  /** If true, perform a grid search for bulk solvent parameters. */
  private boolean gridSearch;
  /** If true, fit a scaling spline between Fo and Fc. */
  private boolean splineFit;

  /**
   * construct a diffraction data assembly, assumes an X-ray data set with a weight of 1.0 using the
   * same name as the molecular assembly
   *
   * @param assembly {@link ffx.potential.MolecularAssembly molecular assembly} object, used as the
   *     atomic model for comparison against the data
   * @param properties system properties file
   */
  public DiffractionData(MolecularAssembly assembly, CompositeConfiguration properties) {
    this(new MolecularAssembly[] {assembly}, properties, POLYNOMIAL, new DiffractionFile(assembly));
  }

  /**
   * construct a diffraction data assembly
   *
   * @param assembly {@link ffx.potential.MolecularAssembly molecular assembly} object, used as the
   *     atomic model for comparison against the data
   * @param properties system properties file
   * @param datafile one or more {@link ffx.xray.parsers.DiffractionFile} to be refined against
   */
  public DiffractionData(
      MolecularAssembly assembly, CompositeConfiguration properties, DiffractionFile... datafile) {
    this(new MolecularAssembly[] {assembly}, properties, POLYNOMIAL, datafile);
  }

  /**
   * construct a diffraction data assembly, assumes an X-ray data set with a weight of 1.0 using the
   * same name as the molecular assembly
   *
   * @param assembly {@link ffx.potential.MolecularAssembly molecular assembly} object, used as the
   *     atomic model for comparison against the data
   * @param properties system properties file
   * @param solventmodel the type of solvent model desired - see {@link
   *     CrystalReciprocalSpace.SolventModel bulk solvent model} selections
   */
  public DiffractionData(
      MolecularAssembly assembly, CompositeConfiguration properties, SolventModel solventmodel) {
    this(
        new MolecularAssembly[] {assembly},
        properties,
        solventmodel,
        new DiffractionFile(assembly));
  }

  /**
   * construct a diffraction data assembly
   *
   * @param assembly {@link ffx.potential.MolecularAssembly molecular assembly} object, used as the
   *     atomic model for comparison against the data
   * @param properties system properties file
   * @param solventmodel the type of solvent model desired - see {@link
   *     CrystalReciprocalSpace.SolventModel bulk solvent model} selections
   * @param datafile one or more {@link ffx.xray.parsers.DiffractionFile} to be refined against
   */
  public DiffractionData(
      MolecularAssembly assembly,
      CompositeConfiguration properties,
      SolventModel solventmodel,
      DiffractionFile... datafile) {
    this(new MolecularAssembly[] {assembly}, properties, solventmodel, datafile);
  }

  /**
   * construct a diffraction data assembly, assumes an X-ray data set with a weight of 1.0 using the
   * same name as the molecular assembly
   *
   * @param assembly {@link ffx.potential.MolecularAssembly molecular assembly} object array
   *     (typically containing alternate conformer assemblies), used as the atomic model for
   *     comparison against the data
   * @param properties system properties file
   */
  public DiffractionData(MolecularAssembly[] assembly, CompositeConfiguration properties) {
    this(assembly, properties, POLYNOMIAL, new DiffractionFile(assembly[0]));
  }

  /**
   * construct a diffraction data assembly
   *
   * @param assembly {@link ffx.potential.MolecularAssembly molecular assembly} object array
   *     (typically containing alternate conformer assemblies), used as the atomic model for
   *     comparison against the data
   * @param properties system properties file
   * @param datafile one or more {@link ffx.xray.parsers.DiffractionFile} to be refined against
   */
  public DiffractionData(
      MolecularAssembly[] assembly,
      CompositeConfiguration properties,
      DiffractionFile... datafile) {
    this(assembly, properties, POLYNOMIAL, datafile);
  }

  /**
   * construct a diffraction data assembly
   *
   * @param assembly {@link ffx.potential.MolecularAssembly molecular assembly} object array
   *     (typically containing alternate conformer assemblies), used as the atomic model for
   *     comparison against the data
   * @param properties system properties file
   * @param solventmodel the type of solvent model desired - see {@link
   *     CrystalReciprocalSpace.SolventModel bulk solvent model} selections
   * @param datafile one or more {@link ffx.xray.parsers.DiffractionFile} to be refined against
   */
  public DiffractionData(
      MolecularAssembly[] assembly,
      CompositeConfiguration properties,
      SolventModel solventmodel,
      DiffractionFile... datafile) {

    this.assembly = assembly;
    this.solventModel = solventmodel;
    this.modelName = assembly[0].getFile().getName();
    this.dataFiles = datafile;
    this.n = datafile.length;

    int rflag = properties.getInt("rfreeflag", -1);
    fsigfCutoff = properties.getDouble("fsigfcutoff", -1.0);
    gridSearch = properties.getBoolean("gridsearch", false);
    splineFit = properties.getBoolean("splinefit", true);
    use_3g = properties.getBoolean("use_3g", true);
    aRadBuff = properties.getDouble("aradbuff", 0.75);
    double sampling = properties.getDouble("sampling", 0.6);
    xrayScaleTol = properties.getDouble("xrayscaletol", 1e-4);
    sigmaATol = properties.getDouble("sigmaatol", 0.05);
    xWeight = properties.getDouble("xweight", 1.0);
    bSimWeight = properties.getDouble("bsimweight", 1.0);
    bNonZeroWeight = properties.getDouble("bnonzeroweight", 1.0);
    bMass = properties.getDouble("bmass", 5.0);
    residueBFactor = properties.getBoolean("residuebfactor", false);
    nResidueBFactor = properties.getInt("nresiduebfactor", 1);
    addAnisou = properties.getBoolean("addanisou", false);
    refineMolOcc = properties.getBoolean("refinemolocc", false);
    occMass = properties.getDouble("occmass", 10.0);

    ForceField forceField = assembly[0].getForceField();
    nativeEnvironmentApproximation =
        forceField.getBoolean("NATIVE_ENVIRONMENT_APPROXIMATION", false);

    crystal = new Crystal[n];
    resolution = new Resolution[n];
    reflectionList = new ReflectionList[n];
    refinementData = new DiffractionRefinementData[n];
    scaleBulkMinimize = new ScaleBulkMinimize[n];
    sigmaAMinimize = new SigmaAMinimize[n];
    splineMinimize = new SplineMinimize[n];
    crystalStats = new CrystalStats[n];

    // Read in Fo/sigFo/FreeR
    File tmp;
    Crystal crystalinit = Crystal.checkProperties(properties);
    Resolution resolutioninit = Resolution.checkProperties(properties);
    if (crystalinit == null || resolutioninit == null) {
      for (int i = 0; i < n; i++) {
        tmp = new File(datafile[i].getFilename());
        reflectionList[i] = datafile[i].getDiffractionfilter().getReflectionList(tmp, properties);

        if (reflectionList[i] == null) {
          logger.info(" Crystal information from the PDB or property file will be used.");
          crystalinit = assembly[i].getCrystal().getUnitCell();
          double res = datafile[i].getDiffractionfilter().getResolution(tmp, crystalinit);
          if (res < 0.0) {
            logger.severe("MTZ/CIF/CNS file does not contain full crystal information!");
          } else {
            resolutioninit = new Resolution(res, sampling);
            reflectionList[i] = new ReflectionList(crystalinit, resolutioninit, properties);
          }
        }
      }
    } else {
      for (int i = 0; i < n; i++) {
        reflectionList[i] = new ReflectionList(crystalinit, resolutioninit, properties);
      }
    }

    for (int i = 0; i < n; i++) {
      crystal[i] = reflectionList[i].crystal;
      resolution[i] = reflectionList[i].resolution;
      refinementData[i] = new DiffractionRefinementData(properties, reflectionList[i]);
      tmp = new File(datafile[i].getFilename());
      datafile[i]
          .getDiffractionfilter()
          .readFile(tmp, reflectionList[i], refinementData[i], properties);
    }

    if (!crystal[0].equals(assembly[0].getCrystal().getUnitCell())) {
      logger.severe(
          "PDB and reflection file crystal information do not match! (check CRYST1 record?)");
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
      sb.append("  Reflections\n");
      sb.append("   F/sigF cutoff: ").append(fsigfCutoff).append("\n");
      sb.append("   R Free flag (-1 auto-determine from the data): ").append(rflag).append("\n");
      sb.append("   Number of bins: ").append(reflectionList[0].nbins).append("\n");
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
    refinementModel = new RefinementModel(assembly, refineMolOcc);

    // initialize atomic form factors
    for (Atom a : refinementModel.getTotalAtomArray()) {
      a.setFormFactorIndex(-1);
      XRayFormFactor atomff = new XRayFormFactor(a, use_3g, 2.0);
      a.setFormFactorIndex(atomff.ffIndex);

      if (a.getOccupancy() == 0.0) {
        a.setFormFactorWidth(1.0);
        continue;
      }

      double arad;
      try {
        arad = a.getVDWType().radius * 0.5;
      } catch (NullPointerException ex) {
        logger.warning(
            format(
                " Failure to get van der Waals type for atom %s; ensure the vdW term is enabled!",
                a.toString()));
        throw ex;
      }
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
    crystalReciprocalSpacesFc = new CrystalReciprocalSpace[n];
    crystalReciprocalSpacesFs = new CrystalReciprocalSpace[n];

    parallelTeam = assembly[0].getParallelTeam();

    String gridString = properties.getString("grid-method", "SLICE").toUpperCase();
    try {
      gridMethod = CrystalReciprocalSpace.GridMethod.valueOf(gridString);
    } catch (Exception e) {
      gridMethod = CrystalReciprocalSpace.GridMethod.SLICE;
    }

    parallelTeam = new ParallelTeam();
    for (int i = 0; i < n; i++) {
      // Atomic Scattering
      crystalReciprocalSpacesFc[i] =
          new CrystalReciprocalSpace(
              reflectionList[i],
              refinementModel.getTotalAtomArray(),
              parallelTeam,
              parallelTeam,
              false,
              dataFiles[i].isNeutron(),
              SolventModel.NONE,
              gridMethod);
      refinementData[i].setCrystalReciprocalSpaceFc(crystalReciprocalSpacesFc[i]);
      crystalReciprocalSpacesFc[i].setUse3G(use_3g);
      crystalReciprocalSpacesFc[i].setWeight(dataFiles[i].getWeight());
      crystalReciprocalSpacesFc[i].lambdaTerm = false;
      crystalReciprocalSpacesFc[i].setNativeEnvironmentApproximation(
          nativeEnvironmentApproximation);

      // Bulk Solvent Scattering
      crystalReciprocalSpacesFs[i] =
          new CrystalReciprocalSpace(
              reflectionList[i],
              refinementModel.getTotalAtomArray(),
              parallelTeam,
              parallelTeam,
              true,
              dataFiles[i].isNeutron(),
              solventmodel,
              gridMethod);
      refinementData[i].setCrystalReciprocalSpaceFs(crystalReciprocalSpacesFs[i]);
      crystalReciprocalSpacesFs[i].setUse3G(use_3g);
      crystalReciprocalSpacesFs[i].setWeight(dataFiles[i].getWeight());
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
   *     (typically containing alternate conformer assemblies), used as the atomic model to average
   *     in with previous assembly
   * @param index the current data index (for cumulative average purposes)
   */
  public void AverageFc(MolecularAssembly[] assembly, int index) {
    RefinementModel tmprefinementmodel = new RefinementModel(assembly, refineMolOcc);

    // initialize atomic form factors
    for (Atom a : tmprefinementmodel.getTotalAtomArray()) {
      a.setFormFactorIndex(-1);
      XRayFormFactor atomff = new XRayFormFactor(a, use_3g, 2.0);
      a.setFormFactorIndex(atomff.ffIndex);

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
              tmprefinementmodel.getTotalAtomArray(),
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
              tmprefinementmodel.getTotalAtomArray(),
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
   * getActiveAtomArray.
   *
   * @return an array of {@link ffx.potential.bonded.Atom} objects.
   */
  public Atom[] getActiveAtomArray() {
    return getAtomArray();
  }

  /** {@inheritDoc} */
  @Override
  public List<List<Molecule>> getAltMolecules() {
    return refinementModel.getAltMolecules();
  }

  /** {@inheritDoc} */
  @Override
  public List<List<Residue>> getAltResidues() {
    return refinementModel.getAltResidues();
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
   * {@inheritDoc}
   *
   * <p>return the atom array for the model associated with this data
   */
  @Override
  public Atom[] getAtomArray() {
    return refinementModel.getTotalAtomArray();
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
   * Getter for the field <code>crystalStats</code>.
   *
   * @return the crystalStats
   */
  public CrystalStats[] getCrystalStats() {
    return crystalStats;
  }

  /**
   * Getter for the field <code>dataFiles</code>.
   *
   * @return the dataFiles
   */
  public DiffractionFile[] getDataFiles() {
    return dataFiles;
  }

  /**
   * Getter for the field <code>fsigfCutoff</code>.
   *
   * @return the fsigfCutoff
   */
  public double getFsigfCutoff() {
    return fsigfCutoff;
  }

  /**
   * Getter for the field <code>modelName</code>.
   *
   * @return the modelName
   */
  public String getModelName() {
    return modelName;
  }

  /** {@inheritDoc} */
  @Override
  public MolecularAssembly[] getMolecularAssemblies() {
    return assembly;
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
   * Return R value for OST x-ray minimization.
   *
   * @return a double.
   */
  public double getRCrystalStat() {
    return crystalStats[0].getR();
  }

  /**
   * Getter for the field <code>refinementData</code>.
   *
   * @return the refinementData
   */
  public DiffractionRefinementData[] getRefinementData() {
    return refinementData;
  }

  /** {@inheritDoc} */
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
   * Getter for the field <code>scaleBulkMinimize</code>.
   *
   * @return the scaleBulkMinimize
   */
  public ScaleBulkMinimize[] getScaleBulkMinimize() {
    return scaleBulkMinimize;
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
   * Getter for the field <code>sigmaAMinimize</code>.
   *
   * @return the sigmaAMinimize
   */
  public SigmaAMinimize[] getSigmaAMinimize() {
    return sigmaAMinimize;
  }

  /**
   * Getter for the field <code>sigmaATol</code>.
   *
   * @return the sigmaATol
   */
  public double getSigmaATol() {
    return sigmaATol;
  }

  /**
   * Getter for the field <code>solventModel</code>.
   *
   * @return the solventModel
   */
  public SolventModel getSolventModel() {
    return solventModel;
  }

  /**
   * Getter for the field <code>splineMinimize</code>.
   *
   * @return the splineMinimize
   */
  public SplineMinimize[] getSplineMinimize() {
    return splineMinimize;
  }

  /** {@inheritDoc} */
  @Override
  public double getWeight() {
    return xWeight;
  }

  /** {@inheritDoc} */
  @Override
  public void setWeight(double weight) {
    this.xWeight = weight;
  }

  /**
   * Getter for the field <code>xrayScaleTol</code>.
   *
   * @return the xrayScaleTol
   */
  public double getXrayScaleTol() {
    return xrayScaleTol;
  }

  /**
   * Getter for the field <code>aRadBuff</code>.
   *
   * @return the aRadBuff
   */
  public double getaRadBuff() {
    return aRadBuff;
  }

  /**
   * Getter for the field <code>xWeight</code>.
   *
   * @return the xWeight
   */
  public double getxWeight() {
    return xWeight;
  }

  /**
   * isRefineMolOcc.
   *
   * @return the refineMolOcc
   */
  public boolean isRefineMolOcc() {
    return refineMolOcc;
  }

  /**
   * isUse_3g.
   *
   * @return the use_3g
   */
  public boolean isUse_3g() {
    return use_3g;
  }

  /** {@inheritDoc} */
  @Override
  public String printEnergyUpdate() {
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < n; i++) {
      sb.append(
          format(
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

  /** {@inheritDoc} */
  @Override
  public String printOptimizationHeader() {
    return "R  Rfree";
  }

  /** {@inheritDoc} */
  @Override
  public String printOptimizationUpdate() {
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < n; i++) {
      sb.append(format("%6.2f %6.2f ", crystalStats[i].getR(), crystalStats[i].getRFree()));
    }

    return sb.toString();
  }

  /** Print scale and R statistics for all datasets associated with the model. */
  public void printScaleAndR() {
    for (int i = 0; i < n; i++) {
      if (!scaled[i]) {
        scaleBulkFit(i);
      }
      crystalStats[i].printScaleStats();
      crystalStats[i].printRStats();
    }
  }

  /** Print all statistics for all datasets associated with the model. */
  public void printStats() {
    int nat = 0;
    int nnonh = 0;

    // Note - we are including inactive and/or un-used atoms in the following loop.
    for (Atom a : refinementModel.getTotalAtomList()) {
      if (a.getOccupancy() == 0.0) {
        continue;
      }
      nat++;
      if (a.getAtomicNumber() == 1) {
        continue;
      }
      nnonh++;
    }

    for (int i = 0; i < n; i++) {
      if (!scaled[i]) {
        scaleBulkFit(i);
      }

      String sb =
          format(
              " Statistics for Data Set %d of %d\n\n"
                  + "  Weight:     %6.2f\n  Neutron data: %4s\n"
                  + "  Model:        %s\n  Data file:    %s\n",
              i + 1,
              n,
              dataFiles[i].getWeight(),
              dataFiles[i].isNeutron(),
              modelName,
              dataFiles[i].getFilename());
      logger.info(sb);

      crystalStats[i].printScaleStats();
      crystalStats[i].printDPIStats(nnonh, nat);
      crystalStats[i].printHKLStats();
      crystalStats[i].printSNStats();
      crystalStats[i].printRStats();
    }
  }

  /** Scale model and fit bulk solvent to all data. */
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
    String sb =
        format(
            " Scaling Data Set %d of %d\n\n  Weight: %6.2f\n  Neutron data: %s\n  Model: %s\n  Data file: %s",
            i + 1,
            n,
            dataFiles[i].getWeight(),
            dataFiles[i].isNeutron(),
            modelName,
            dataFiles[i].getFilename());
    logger.info(sb);

    // reset some values
    refinementData[i].bulkSolventK = 0.33;
    refinementData[i].bulkSolventUeq = 50.0 / (8.0 * Math.PI * Math.PI);
    refinementData[i].modelScaleK = 0.0;
    for (int j = 0; j < 6; j++) {
      refinementData[i].modelAnisoB[j] = 0.0;
    }

    // run FFTs
    crystalReciprocalSpacesFc[i].computeDensity(refinementData[i].fc);
    if (solventModel != SolventModel.NONE) {
      crystalReciprocalSpacesFs[i].computeDensity(refinementData[i].fs);
    }

    // initialize minimizers
    scaleBulkMinimize[i] =
        new ScaleBulkMinimize(
            reflectionList[i], refinementData[i], crystalReciprocalSpacesFs[i], parallelTeam);
    splineMinimize[i] =
        new SplineMinimize(
            reflectionList[i], refinementData[i], refinementData[i].spline, SplineEnergy.Type.FOFC);

    // minimize
    if (solventModel != SolventModel.NONE && gridSearch) {
      scaleBulkMinimize[i].minimize(6, xrayScaleTol);
      scaleBulkMinimize[i].gridOptimize();
    }
    scaleBulkMinimize[i].minimize(6, xrayScaleTol);

    // sigmaA / LLK calculation
    sigmaAMinimize[i] = new SigmaAMinimize(reflectionList[i], refinementData[i], parallelTeam);
    sigmaAMinimize[i].minimize(7, sigmaATol);

    if (splineFit) {
      splineMinimize[i].minimize(7, 1e-5);
    }

    scaled[i] = true;
  }

  /**
   * set the bulk solvent parameters for a given bulk solvent model
   *
   * @param a typically the width of the atom
   * @param b typically the rate with which the atom transitions to bulk
   * @see CrystalReciprocalSpace#setSolventAB(double, double)
   */
  public void setSolventAB(double a, double b) {
    for (int i = 0; i < n; i++) {
      if (solventModel != SolventModel.NONE) {
        refinementData[i].solventA = a;
        refinementData[i].solventB = b;
        crystalReciprocalSpacesFs[i].setSolventAB(a, b);
      }
    }
  }

  /** Perform 10 Fc calculations for the purposes of timings. */
  public void timings() {
    logger.info(" Performing 10 Fc calculations for timing...");
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < 10; j++) {
        crystalReciprocalSpacesFc[i].computeDensity(refinementData[i].fc, true);
        crystalReciprocalSpacesFs[i].computeDensity(refinementData[i].fs, true);
        crystalReciprocalSpacesFc[i].computeAtomicGradients(
            refinementData[i].dFc,
            refinementData[i].freeR,
            refinementData[i].rFreeFlag,
            RefinementMode.COORDINATES,
            true);
        crystalReciprocalSpacesFs[i].computeAtomicGradients(
            refinementData[i].dFs,
            refinementData[i].freeR,
            refinementData[i].rFreeFlag,
            RefinementMode.COORDINATES,
            true);
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
        writeData("" + FilenameUtils.removeExtension(filename) + "_" + i + ".mtz", i);
      }
    }
  }

  /**
   * write dataset i to MTZ file
   *
   * @param filename output filename
   * @param i dataset to write out
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
  public void writeMaps(String filename) {
    if (n == 1) {
      writeMaps(filename, 0);
    } else {
      for (int i = 0; i < n; i++) {
        writeMaps("" + FilenameUtils.removeExtension(filename) + "_" + i + ".map", i);
      }
    }
  }

  /**
   * Write current model to PDB file.
   *
   * @param filename output PDB filename
   */
  public void writeModel(String filename) {
    StringBuilder remark = new StringBuilder();
    File file = new File(filename);
    PDBFilter pdbFilter = new PDBFilter(file, Arrays.asList(assembly), null, null);

    Date now = new Date();
    SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss ");
    remark.append("REMARK FFX output ISO-8601 date: " + sdf.format(now) + "\n");
    remark.append("REMARK\n");
    remark.append("REMARK   3\n");
    remark.append("REMARK   3 REFINEMENT\n");
    remark.append("REMARK   3   PROGRAM     : FORCE FIELD X\n");
    remark.append("REMARK   3\n");
    for (int i = 0; i < n; i++) {
      remark.append("REMARK   3  DATA SET " + (i + 1) + "\n");
      if (dataFiles[i].isNeutron()) {
        remark.append("REMARK   3   DATA SET TYPE   : NEUTRON\n");
      } else {
        remark.append("REMARK   3   DATA SET TYPE   : X-RAY\n");
      }
      remark.append("REMARK   3   DATA SET WEIGHT : " + dataFiles[i].getWeight() + "\n");
      remark.append("REMARK   3\n");
      remark.append(crystalStats[i].getPDBHeaderString());
    }
    for (int i = 0; i < assembly.length; i++) {
      remark.append("REMARK   3  CHEMICAL SYSTEM " + (i + 1) + "\n");
      remark.append(assembly[i].getPotentialEnergy().getPDBHeaderString());
    }
    pdbFilter.writeFileWithHeader(file, remark);
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
        writeSolventMask("" + FilenameUtils.removeExtension(filename) + "_" + i + ".map", i);
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
        writeSolventMaskCNS("" + FilenameUtils.removeExtension(filename) + "_" + i + ".map", i);
      }
    }
  }

  /**
   * Move the atomic coordinates for the FFT calculation - this is independent of the {@link
   * ffx.potential.bonded.Atom} object coordinates as it uses a linearized array
   *
   * @param x array of coordinates to move atoms to
   * @see CrystalReciprocalSpace#setCoordinates(double[])
   */
  void setFFTCoordinates(double[] x) {
    for (int i = 0; i < n; i++) {
      crystalReciprocalSpacesFc[i].setCoordinates(x);
      crystalReciprocalSpacesFs[i].setCoordinates(x);
    }
  }

  /**
   * Determine the total atomic gradients due to the diffraction data - performs an inverse FFT
   *
   * @param refinementMode the {@link RefinementMinimize.RefinementMode refinement mode} requested
   * @see CrystalReciprocalSpace#computeAtomicGradients(double[][], int[], int,
   *     ffx.xray.RefinementMinimize.RefinementMode, boolean) computeAtomicGradients
   */
  void computeAtomicGradients(RefinementMode refinementMode) {
    for (int i = 0; i < n; i++) {
      crystalReciprocalSpacesFc[i].computeAtomicGradients(
          refinementData[i].dFc,
          refinementData[i].freeR,
          refinementData[i].rFreeFlag,
          refinementMode);
      crystalReciprocalSpacesFs[i].computeAtomicGradients(
          refinementData[i].dFs,
          refinementData[i].freeR,
          refinementData[i].rFreeFlag,
          refinementMode);
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
   * write 2Fo-Fc and Fo-Fc maps for a datasets
   *
   * @param filename output root filename for Fo-Fc and 2Fo-Fc maps
   * @param i a int.
   */
  private void writeMaps(String filename, int i) {
    if (!scaled[i]) {
      scaleBulkFit(i);
    }

    // Fo-Fc
    crystalReciprocalSpacesFc[i].computeAtomicGradients(
        refinementData[i].foFc1,
        refinementData[i].freeR,
        refinementData[i].rFreeFlag,
        RefinementMode.COORDINATES);
    double[] densityGrid = crystalReciprocalSpacesFc[i].getDensityGrid();
    int extx = (int) crystalReciprocalSpacesFc[i].getXDim();
    int exty = (int) crystalReciprocalSpacesFc[i].getYDim();
    int extz = (int) crystalReciprocalSpacesFc[i].getZDim();

    CCP4MapWriter mapwriter =
        new CCP4MapWriter(
            extx, exty, extz, crystal[i], FilenameUtils.removeExtension(filename) + "_fofc.map");
    mapwriter.write(densityGrid);

    // 2Fo-Fc
    crystalReciprocalSpacesFc[i].computeAtomicGradients(
        refinementData[i].foFc2,
        refinementData[i].freeR,
        refinementData[i].rFreeFlag,
        RefinementMode.COORDINATES);
    densityGrid = crystalReciprocalSpacesFc[i].getDensityGrid();
    extx = (int) crystalReciprocalSpacesFc[i].getXDim();
    exty = (int) crystalReciprocalSpacesFc[i].getYDim();
    extz = (int) crystalReciprocalSpacesFc[i].getZDim();

    mapwriter =
        new CCP4MapWriter(
            extx, exty, extz, crystal[i], FilenameUtils.removeExtension(filename) + "_2fofc.map");
    mapwriter.write(densityGrid);
  }

  /**
   * Write bulk solvent mask for dataset i to a CNS map file.
   *
   * @param filename output filename
   * @param i dataset to write out
   */
  private void writeSolventMaskCNS(String filename, int i) {
    if (solventModel != SolventModel.NONE) {
      try {
        PrintWriter cnsfile = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
        cnsfile.println(" ANOMalous=FALSE");
        cnsfile.println(" DECLare NAME=FS DOMAin=RECIprocal TYPE=COMP END");
        for (HKL ih : reflectionList[i].hkllist) {
          int j = ih.index();
          cnsfile.printf(
              " INDE %d %d %d FS= %.4f %.4f\n",
              ih.h(),
              ih.k(),
              ih.l(),
              refinementData[i].fsF(j),
              Math.toDegrees(refinementData[i].fsPhi(j)));
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
   * @param i dataset to write out
   * @see CCP4MapWriter#write(double[], boolean)
   */
  private void writeSolventMask(String filename, int i) {
    if (solventModel != SolventModel.NONE) {
      CCP4MapWriter mapwriter =
          new CCP4MapWriter(
              (int) crystalReciprocalSpacesFs[i].getXDim(),
              (int) crystalReciprocalSpacesFs[i].getYDim(),
              (int) crystalReciprocalSpacesFs[i].getZDim(),
              crystal[i],
              filename);
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

  /**
   * Getter for the field <code>bMass</code>.
   *
   * @return the bMass
   */
  double getbMass() {
    return bMass;
  }

  /**
   * isResidueBFactor.
   *
   * @return the residueBFactor
   */
  boolean isResidueBFactor() {
    return residueBFactor;
  }

  /**
   * Getter for the field <code>nResidueBFactor</code>.
   *
   * @return the nResidueBFactor
   */
  int getnResidueBFactor() {
    return nResidueBFactor;
  }

  /**
   * isAddAnisou.
   *
   * @return the addAnisou
   */
  boolean isAddAnisou() {
    return addAnisou;
  }

  /**
   * Getter for the field <code>occMass</code>.
   *
   * @return the occMass
   */
  double getOccMass() {
    return occMass;
  }
}
