/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.xray;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.apache.commons.configuration.CompositeConfiguration;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.Keyword;

/**
 * @author Timothy D. Fenn
 */
@RunWith(Parameterized.class)
public class ScaleBulkMinimizeTest {

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {false,
                "NSF D2 domain test",
                "ffx/xray/structures/1NSF.pdb",
                "ffx/xray/structures/1NSF.mtz",
                null,
                25.178605682089,
                25.448314839595,
                0.8921390108139,
                0.1526816244814},
            {false,
                "SNARE complex",
                "ffx/xray/structures/1N7S.pdb",
                "ffx/xray/structures/1N7S.mtz",
                null,
                19.412671496011,
                21.555930987573,
                0.9314271250347,
                0.1361239311856}
        });
    }
    private final String info;
    private final CrystalStats crystalStats;
    private final double r;
    private final double rFree;
    private final double sigmaA;
    private final double sigmaW;
    private final boolean ci;
    private final boolean ciOnly;

    public ScaleBulkMinimizeTest(boolean ciOnly,
            String info, String pdbname, String mtzname, String cifname,
            double r, double rFree, double sigmaA, double sigmaW) {
        this.ciOnly = ciOnly;
        this.info = info;
        this.r = r;
        this.rFree = rFree;
        this.sigmaA = sigmaA;
        this.sigmaW = sigmaW;

        String ffxCi = System.getProperty("ffx.ci");
        if (ffxCi != null && ffxCi.equalsIgnoreCase("true")) {
            ci = true;
        } else {
            ci = false;
        }

        if (!ci && ciOnly) {
            crystalStats = null;
            return;
        }

        int index = pdbname.lastIndexOf(".");
        String name = pdbname.substring(0, index);

        // load the structure
        ClassLoader cl = this.getClass().getClassLoader();
        File structure = new File(cl.getResource(pdbname).getPath());
        File mtzFile = null, cifFile = null;
        if (mtzname != null) {
            mtzFile = new File(cl.getResource(mtzname).getPath());
        } else {
            cifFile = new File(cl.getResource(cifname).getPath());
        }

        // load any properties associated with it
        CompositeConfiguration properties = Keyword.loadProperties(structure);

        // read in Fo/sigFo/FreeR
        MTZFilter mtzFilter = new MTZFilter();
        CIFFilter cifFilter = new CIFFilter();
        ReflectionList reflectionList;
        Crystal crystal = Crystal.checkProperties(properties);
        Resolution resolution = Resolution.checkProperties(properties);
        if (crystal == null || resolution == null) {
            if (mtzname != null) {
                reflectionList = mtzFilter.getReflectionList(mtzFile);
            } else {
                reflectionList = cifFilter.getReflectionList(cifFile);
            }
        } else {
            reflectionList = new ReflectionList(crystal, resolution);
        }

        DiffractionRefinementData refinementData = new DiffractionRefinementData(properties,
                reflectionList);
        if (mtzname != null) {
            assertTrue(info + " mtz file should be read in without errors",
                    mtzFilter.readFile(mtzFile, reflectionList, refinementData,
                            properties));
        } else {
            assertTrue(info + " cif file should be read in without errors",
                    cifFilter.readFile(cifFile, reflectionList, refinementData,
                            properties));
        }

        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
        ForceField forceField = forceFieldFilter.parse();

        // associate molecular assembly with the structure, set up forcefield
        MolecularAssembly molecularAssembly = new MolecularAssembly(name);
        molecularAssembly.setFile(structure);
        molecularAssembly.setForceField(forceField);
        PDBFilter pdbFile = new PDBFilter(structure, molecularAssembly, forceField, properties);
        pdbFile.readFile();
        molecularAssembly.finalize(true, forceField);
        ForceFieldEnergy energy = new ForceFieldEnergy(molecularAssembly);

        List<Atom> atomList = molecularAssembly.getAtomList();
        Atom atomArray[] = atomList.toArray(new Atom[atomList.size()]);

        // set up FFT and run it
        ParallelTeam parallelTeam = new ParallelTeam();
        CrystalReciprocalSpace crs = new CrystalReciprocalSpace(reflectionList,
                atomArray, parallelTeam, parallelTeam, false);
        crs.computeDensity(refinementData.fc);
        refinementData.setCrystalReciprocalSpace_fc(crs);
        crs = new CrystalReciprocalSpace(reflectionList,
                atomArray, parallelTeam, parallelTeam, true);
        crs.computeDensity(refinementData.fs);
        refinementData.setCrystalReciprocalSpace_fs(crs);

        ScaleBulkMinimize scaleBulkMinimize
                = new ScaleBulkMinimize(reflectionList, refinementData, crs, parallelTeam);
        scaleBulkMinimize.minimize(6, 1.0e-4);

        SigmaAMinimize sigmaAMinimize = new SigmaAMinimize(reflectionList,
                refinementData, parallelTeam);
        sigmaAMinimize.minimize(7, 1.0e-1);

        SplineMinimize splineMinimize = new SplineMinimize(reflectionList,
                refinementData, refinementData.spline, SplineEnergy.Type.FOFC);
        splineMinimize.minimize(7, 1e-5);

        crystalStats = new CrystalStats(reflectionList, refinementData);
    }

    @Test
    public void testCrystalStats() {
        if (!ci && ciOnly) {
            return;
        }

        crystalStats.printScaleStats();
        crystalStats.printHKLStats();
        crystalStats.printSNStats();
        crystalStats.printRStats();

        assertEquals(info + " R value",
                r, crystalStats.getR(), 0.01);
        assertEquals(info + " Rfree value",
                rFree, crystalStats.getRFree(), 0.01);
        assertEquals(info + " sigmaA s",
                sigmaA, crystalStats.getSigmaA(), 0.001);
        assertEquals(info + " sigmaA w",
                sigmaW, crystalStats.getSigmaW(), 0.001);
    }
}
