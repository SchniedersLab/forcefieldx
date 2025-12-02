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

import static java.lang.String.format;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import edu.rit.pj.ParallelTeam;
import ffx.algorithms.misc.AlgorithmsTest;
import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.utils.PotentialsUtils;
import ffx.xray.parsers.CIFFilter;
import ffx.xray.parsers.MTZFilter;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.logging.Logger;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/**
 * @author Timothy D. Fenn and Michael J. Schnieders
 */
@RunWith(Parameterized.class)
public class XRayMinimizeTest extends AlgorithmsTest {

  private static final Logger logger = Logger.getLogger(XRayMinimizeTest.class.getName());
  private final String info;
  private final double r;
  private final double rFree;
  private final double sigmaA;
  private final double sigmaW;
  private final boolean ciOnly;
  private DiffractionRefinementData refinementData;
  private ReflectionList reflectionList;
  private ParallelTeam parallelTeam;

  public XRayMinimizeTest(boolean ciOnly, String info, String pdbname, String mtzname, String cifname,
                          double r, double rFree, double sigmaA, double sigmaW) {
    this.ciOnly = ciOnly;
    this.info = info;
    this.r = r;
    this.rFree = rFree;
    this.sigmaA = sigmaA;
    this.sigmaW = sigmaW;

    if (!ffxCI && ciOnly) {
      return;
    }

    // Load the structure
    File structure = getResourceFile(pdbname);
    File mtzFile = getResourceFile(mtzname);
    File cifFile = getResourceFile(cifname);
    PotentialsUtils potutil = new PotentialsUtils();
    MolecularAssembly mola = potutil.open(structure);
    CompositeConfiguration properties = mola.getProperties();

    // read in Fo/sigFo/FreeR
    MTZFilter mtzFilter = new MTZFilter();
    CIFFilter cifFilter = new CIFFilter();
    Crystal crystal = Crystal.checkProperties(properties);
    double defaultResolution = -1.0;
    if (crystal != null) {
      defaultResolution = mtzFilter.getResolution(mtzFile, crystal);
    }
    Resolution resolution = Resolution.checkProperties(properties, false, defaultResolution);
    if (crystal == null || resolution == null) {
      if (mtzname != null) {
        reflectionList = mtzFilter.getReflectionList(mtzFile);
      } else {
        reflectionList = cifFilter.getReflectionList(cifFile);
      }
    } else {
      reflectionList = new ReflectionList(crystal, resolution);
    }

    refinementData = new DiffractionRefinementData(properties, reflectionList);
    if (mtzname != null) {
      assertTrue(info + " mtz file should be read in without errors",
          mtzFilter.readFile(mtzFile, reflectionList, refinementData, properties));
    } else {
      assertTrue(info + " cif file should be read in without errors",
          cifFilter.readFile(cifFile, reflectionList, refinementData, properties));
    }

    mola.finalize(true, mola.getForceField());

    List<Atom> atomList = mola.getAtomList();
    Atom[] atomArray = atomList.toArray(new Atom[0]);

    // Set up the CrystalReciprocalSpace
    parallelTeam = new ParallelTeam();
    CrystalReciprocalSpace crs =
        new CrystalReciprocalSpace(reflectionList, atomArray, parallelTeam, parallelTeam, false);
    crs.computeDensity(refinementData.fc);
    refinementData.setCrystalReciprocalSpaceFc(crs);
    crs = new CrystalReciprocalSpace(reflectionList, atomArray, parallelTeam, parallelTeam, true);
    crs.computeDensity(refinementData.fs);
    refinementData.setCrystalReciprocalSpaceFs(crs);

    ScaleBulkMinimize scaleBulkMinimize = new ScaleBulkMinimize(reflectionList,
        refinementData, refinementData.crystalReciprocalSpaceFs, parallelTeam);
    scaleBulkMinimize.minimize(6, 1.0e-4);

    SigmaAMinimize sigmaAMinimize = new SigmaAMinimize(reflectionList, refinementData, parallelTeam);
    sigmaAMinimize.minimize(7, 2.0e-2);

    SplineMinimize splineMinimize = new SplineMinimize(reflectionList, refinementData,
        refinementData.spline, SplineEnergy.SplineType.FOFC);
    splineMinimize.minimize(7, 1e-5);
  }

  @Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][]{
            {
                false,
                "NSF D2 domain test",
                "1NSF.pdb",
                "1NSF.mtz",
                null,
                25.073186034370938,
                25.310851476829942,
                0.895228035088474,
                0.1465683355168056
            },
            {
                false,
                "SNARE complex",
                "1N7S.pdb",
                "1N7S.mtz",
                null,
                19.43334760963002,
                21.555930987392596,
                0.9336845537932159,
                0.1319269157669047
            }
        });
  }

  @Test
  public void testCrystalStats() {
    if (!ffxCI && ciOnly) {
      return;
    }
    CrystalStats crystalStats = new CrystalStats(reflectionList, refinementData);

    crystalStats.printScaleStats();
    crystalStats.printHKLStats();
    crystalStats.printSNStats();
    crystalStats.printRStats();

    assertEquals(info + " R value", r, crystalStats.getR(), 0.02);
    assertEquals(info + " Rfree value", rFree, crystalStats.getRFree(), 0.02);
    assertEquals(info + " sigmaA s", sigmaA, crystalStats.getSigmaA(), 0.001);
    assertEquals(info + " sigmaA w", sigmaW, crystalStats.getSigmaW(), 0.001);
  }

  @Test
  public void testScaleBulk() {
    if (!ffxCI && ciOnly) {
      return;
    }
    ScaleBulkMinimize scaleBulkMinimize = new ScaleBulkMinimize(
        reflectionList, refinementData, refinementData.crystalReciprocalSpaceFs, parallelTeam);
    ScaleBulkEnergy scaleBulkEnergy = scaleBulkMinimize.getScaleBulkEnergy();
    int n = scaleBulkMinimize.getNumberOfVariables();
    double[] x = new double[n];
    double[] g = new double[n];
    scaleBulkMinimize.getCoordinates(x);
    scaleBulkEnergy.energyAndGradient(x, g);
    double delta = 1.0e-4;
    double tolerance = 1.0e-4;
    for (int i = 0; i < n; i++) {
      String test = format(" Scale Bulk Solvent Derivative %d.", i);
      double orig = x[i];
      x[i] += delta;
      double ePlus = scaleBulkEnergy.energy(x);
      x[i] -= 2.0 * delta;
      double eMinus = scaleBulkEnergy.energy(x);
      x[i] = orig;
      double fd = (ePlus - eMinus) / (2.0 * delta);
      logger.info(format(" %s A %16.8f vs. FD %16.8f", test, g[i], fd));
      assertEquals(test, g[i], fd, tolerance);
    }
  }

  @Test
  public void testSigmaA() {
    if (!ffxCI && ciOnly) {
      return;
    }
    SigmaAMinimize sigmaAMinimize = new SigmaAMinimize(reflectionList, refinementData, parallelTeam);
    SigmaAEnergy sigmaAEnergy = sigmaAMinimize.getSigmaAEnergy();
    int n = sigmaAMinimize.getNumberOfVariables();
    double[] x = new double[n];
    double[] g = new double[n];
    sigmaAMinimize.getCoordinates(x);
    sigmaAEnergy.energyAndGradient(x, g);
    double delta = 1.0e-4;
    double tolerance = 1.0e-3;
    for (int i = 0; i < n; i++) {
      double orig = x[i];
      x[i] += delta;
      double ePlus = sigmaAEnergy.energy(x);
      x[i] -= 2.0 * delta;
      double eMinus = sigmaAEnergy.energy(x);
      x[i] = orig;
      double fd = (ePlus - eMinus) / (2.0 * delta);
      String test = format(" SigmaA Derivative %d.", i);
      logger.info(format(" %s A %16.8f vs. FD %16.8f", test, g[i], fd));
      assertEquals(test, 1.0, g[i] / fd, tolerance);
    }
  }

  @Test
  public void testSplineFOFC() {
    if (!ffxCI && ciOnly) {
      return;
    }
    SplineMinimize splineMinimize = new SplineMinimize(reflectionList, refinementData,
        refinementData.spline, SplineEnergy.SplineType.FOFC);
    SplineEnergy splineEnergy = splineMinimize.getSplineEnergy();
    int n = splineMinimize.getNumberOfVariables();
    double[] x = new double[n];
    double[] g = new double[n];
    splineMinimize.getCoordinates(x);
    splineEnergy.energyAndGradient(x, g);
    double delta = 1.0e-4;
    double tolerance = 1.0e-4;
    for (int i = 0; i < n; i++) {
      double orig = x[i];
      x[i] += delta;
      double ePlus = splineEnergy.energy(x);
      x[i] -= 2.0 * delta;
      double eMinus = splineEnergy.energy(x);
      x[i] = orig;
      double fd = (ePlus - eMinus) / (2.0 * delta);
      String test = format(" FOFC Spline Derivative %d.", i);
      logger.info(format(" %s A %16.8f vs. FD %16.8f", test, g[i], fd));
      assertEquals(test, 1.0, g[i] / fd, tolerance);
    }
  }

  @Test
  public void testSplineF1F2() {
    if (!ffxCI && ciOnly) {
      return;
    }
    SplineMinimize splineMinimize = new SplineMinimize(
        reflectionList, refinementData, refinementData.spline, SplineEnergy.SplineType.F1F2);
    SplineEnergy splineEnergy = splineMinimize.getSplineEnergy();
    int n = splineMinimize.getNumberOfVariables();
    double[] x = new double[n];
    double[] g = new double[n];
    splineMinimize.getCoordinates(x);
    splineEnergy.energyAndGradient(x, g);
    double delta = 1.0e-4;
    double tolerance = 1.0e-4;
    for (int i = 0; i < n; i++) {
      double orig = x[i];
      x[i] += delta;
      double ePlus = splineEnergy.energy(x);
      x[i] -= 2.0 * delta;
      double eMinus = splineEnergy.energy(x);
      x[i] = orig;
      double fd = (ePlus - eMinus) / (2.0 * delta);
      String test = format(" F1F2 Spline Derivative %d.", i);
      logger.info(format(" %s A %16.8f vs. FD %16.8f", test, g[i], fd));
      assertEquals(test, 1.0, g[i] / fd, tolerance);
    }
  }

  @Test
  public void testSplineFCTOESQ() {
    if (!ffxCI && ciOnly) {
      return;
    }
    SplineMinimize splineMinimize = new SplineMinimize(
        reflectionList, refinementData, refinementData.spline, SplineEnergy.SplineType.FCTOESQ);
    SplineEnergy splineEnergy = splineMinimize.getSplineEnergy();
    int n = splineMinimize.getNumberOfVariables();
    double[] x = new double[n];
    double[] g = new double[n];
    splineMinimize.getCoordinates(x);
    splineEnergy.energyAndGradient(x, g);
    double delta = 1.0e-4;
    double tolerance = 1.0e-4;
    for (int i = 0; i < n; i++) {
      double orig = x[i];
      x[i] += delta;
      double ePlus = splineEnergy.energy(x);
      x[i] -= 2.0 * delta;
      double eMinus = splineEnergy.energy(x);
      x[i] = orig;
      double fd = (ePlus - eMinus) / (2.0 * delta);
      String test = format(" FCTOESQ Spline Derivative %d.", i);
      logger.info(format(" %s A %16.8f vs. FD %16.8f", test, g[i], fd));
      assertEquals(test, 1.0, g[i] / fd, tolerance);
    }
  }

  @Test
  public void testSplineFOTOESQ() {
    if (!ffxCI && ciOnly) {
      return;
    }
    SplineMinimize splineMinimize = new SplineMinimize(
        reflectionList, refinementData, refinementData.spline, SplineEnergy.SplineType.FOTOESQ);
    SplineEnergy splineEnergy = splineMinimize.getSplineEnergy();
    int n = splineMinimize.getNumberOfVariables();
    double[] x = new double[n];
    double[] g = new double[n];
    splineMinimize.getCoordinates(x);
    splineEnergy.energyAndGradient(x, g);
    double delta = 1.0e-4;
    double tolerance = 1.0e-4;
    for (int i = 0; i < n; i++) {
      double orig = x[i];
      x[i] += delta;
      double ePlus = splineEnergy.energy(x);
      x[i] -= 2.0 * delta;
      double eMinus = splineEnergy.energy(x);
      x[i] = orig;
      double fd = (ePlus - eMinus) / (2.0 * delta);
      String test = format(" FOTOESQ Spline Derivative %d.", i);
      logger.info(format(" %s A %16.8f vs. FD %16.8f", test, g[i], fd));
      assertEquals(test, 1.0, g[i] / fd, tolerance);
    }
  }
}
