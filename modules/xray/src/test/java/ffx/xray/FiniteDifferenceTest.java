/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.Keyword;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.RefinementMinimize.RefinementMode;
import ffx.xray.parsers.MTZFilter;

import static ffx.numerics.VectorMath.b2u;
import static ffx.xray.CrystalReciprocalSpace.SolventModel.NONE;

/**
 * @author Timothy D. Fenn
 */
@RunWith(Parameterized.class)
public class FiniteDifferenceTest {

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {false,
                "ala met anisou",
                NONE,
                new int[]{91, 105, 119},
                "ffx/xray/structures/alamet.pdb",
                "ffx/xray/structures/alamet.mtz"}
        });
    }

    private final boolean ci;
    private final boolean ciOnly;
    private final Atom atomArray[];
    private final int atoms[];
    private final DiffractionRefinementData refinementData;
    private final SigmaAMinimize sigmaAMinimize;

    public FiniteDifferenceTest(boolean ciOnly,
            String info, SolventModel solventModel, int[] atoms,
            String pdbName, String mtzName) {
        this.ciOnly = ciOnly;
        this.atoms = atoms;

        ci = System.getProperty("ffx.ci","false").equalsIgnoreCase("true");
        if (!ci && ciOnly) {
            atomArray = null;
            refinementData = null;
            sigmaAMinimize = null;
            return;
        }

        int index = pdbName.lastIndexOf(".");
        String name = pdbName.substring(0, index);

        // load the structure
        ClassLoader cl = this.getClass().getClassLoader();
        File structure = new File(cl.getResource(pdbName).getPath());
        File mtzFile = new File(cl.getResource(mtzName).getPath());

        // load any properties associated with it
        CompositeConfiguration properties = Keyword.loadProperties(structure);

        // read in Fo/sigFo/FreeR
        MTZFilter mtzFilter = new MTZFilter();
        ReflectionList reflectionList;
        Crystal crystal = Crystal.checkProperties(properties);
        Resolution resolution = Resolution.checkProperties(properties);
        if (crystal == null || resolution == null) {
            reflectionList = mtzFilter.getReflectionList(mtzFile);
        } else {
            reflectionList = new ReflectionList(crystal, resolution);
        }

        refinementData = new DiffractionRefinementData(properties,
                reflectionList);
        assertTrue(info + " mtz file should be read in without errors",
                mtzFilter.readFile(mtzFile, reflectionList, refinementData,
                        properties));

        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
        ForceField forceField = forceFieldFilter.parse();

        // associate molecular assembly with the structure, set up forcefield
        MolecularAssembly molecularAssembly = new MolecularAssembly(name);
        molecularAssembly.setFile(structure);
        molecularAssembly.setForceField(forceField);
        PDBFilter pdbFile = new PDBFilter(structure, molecularAssembly, forceField, properties);
        pdbFile.readFile();
        pdbFile.applyAtomProperties();
        molecularAssembly.finalize(true, forceField);
        ForceFieldEnergy energy = new ForceFieldEnergy(molecularAssembly, pdbFile.getCoordRestraints());

        List<Atom> atomList = molecularAssembly.getAtomList();
        atomArray = atomList.toArray(new Atom[0]);
        boolean use_3g = properties.getBoolean("use_3g", true);

        // initialize atomic form factors
        for (Atom atom : atomArray) {
            XRayFormFactor xrayFormFactor = new XRayFormFactor(atom, use_3g, 2.0);
            atom.setFormFactorIndex(xrayFormFactor.ffIndex);
            if (atom.getOccupancy() == 0.0) {
                atom.setFormFactorWidth(1.0);
                continue;
            }
            double aRad = 2.4;
            double xyz[] = new double[3];
            xyz[0] = atom.getX() + aRad;
            xyz[1] = atom.getY();
            xyz[2] = atom.getZ();
            while (true) {
                double rho = xrayFormFactor.rho(0.0, 1.0, xyz);
                if (rho > 0.1) {
                    aRad += 0.5;
                } else if (rho > 0.001) {
                    aRad += 0.1;
                } else {
                    aRad += 0.75;
                    atom.setFormFactorWidth(aRad);
                    break;
                }
                xyz[0] = atom.getX() + aRad;
            }
        }

        // set up FFT and run it
        ParallelTeam parallelTeam = new ParallelTeam();
        CrystalReciprocalSpace crs = new CrystalReciprocalSpace(reflectionList,
                atomArray, parallelTeam, parallelTeam, false, false);
        crs.computeDensity(refinementData.fc);
        refinementData.setCrystalReciprocalSpace_fc(crs);
        crs = new CrystalReciprocalSpace(reflectionList,
                atomArray, parallelTeam, parallelTeam, true, false, solventModel);
        crs.computeDensity(refinementData.fs);
        refinementData.setCrystalReciprocalSpace_fs(crs);

        ScaleBulkMinimize scaleBulkMinimize
                = new ScaleBulkMinimize(reflectionList, refinementData, crs, parallelTeam);
        scaleBulkMinimize.minimize(6, 1.0e-4);

        sigmaAMinimize = new SigmaAMinimize(reflectionList,
                refinementData, parallelTeam);
        sigmaAMinimize.minimize(7, 2.0e-2);

        SplineMinimize splineMinimize = new SplineMinimize(reflectionList,
                refinementData, refinementData.spline,
                SplineEnergy.Type.FOFC);
        splineMinimize.minimize(7, 1e-5);

        CrystalStats crystalstats = new CrystalStats(reflectionList, refinementData);
        crystalstats.printScaleStats();
        crystalstats.printHKLStats();
        crystalstats.printSNStats();
        crystalstats.printRStats();
    }

    @Test
    public void testFiniteDifferences() {
        if (!ci && ciOnly) {
            return;
        }

        // delta for numerical finite differences
        double delta = 1e-4;

        // compute gradients
        refinementData.crs_fc.computeAtomicGradients(refinementData.dfc,
                refinementData.freer, refinementData.rfreeflag,
                RefinementMode.COORDINATES_AND_BFACTORS);

        double mean = 0.0;
        double nmean = 0.0;
        double gxyz[] = new double[3];
        for (int i = 0; i < atoms.length; i++) {
            Atom atom = atomArray[atoms[i]];
            int index = atom.getXYZIndex() - 1;
            if (atom.getOccupancy() == 0.0) {
                continue;
            }

            System.out.println(" " + i + ": " + atom.toString());
            atom.getXYZGradient(gxyz);
            double bg = atom.getTempFactorGradient();

            refinementData.crs_fc.deltaX(index, delta);
            refinementData.crs_fc.computeDensity(refinementData.fc);
            double llk1 = sigmaAMinimize.calculateLikelihood();
            refinementData.crs_fc.deltaX(index, -delta);
            refinementData.crs_fc.computeDensity(refinementData.fc);
            double llk2 = sigmaAMinimize.calculateLikelihood();
            double fd = (llk1 - llk2) / (2.0 * delta);
            System.out.print(String.format(" X A: %16.8f FD: %16.8f Error: %16.8f\n", gxyz[0], fd, gxyz[0] - fd));
            refinementData.crs_fc.deltaX(index, 0.0);

            nmean++;
            mean += (gxyz[0] / fd - mean) / nmean;

            refinementData.crs_fc.deltaY(index, delta);
            refinementData.crs_fc.computeDensity(refinementData.fc);
            llk1 = sigmaAMinimize.calculateLikelihood();
            refinementData.crs_fc.deltaY(index, -delta);
            refinementData.crs_fc.computeDensity(refinementData.fc);
            llk2 = sigmaAMinimize.calculateLikelihood();
            fd = (llk1 - llk2) / (2.0 * delta);
            System.out.print(String.format(" Y A: %16.8f FD: %16.8f Error: %16.8f\n", gxyz[1], fd, gxyz[1] - fd));
            refinementData.crs_fc.deltaY(index, 0.0);

            nmean++;
            mean += (gxyz[1] / fd - mean) / nmean;

            refinementData.crs_fc.deltaZ(index, delta);
            refinementData.crs_fc.computeDensity(refinementData.fc);
            llk1 = sigmaAMinimize.calculateLikelihood();
            refinementData.crs_fc.deltaZ(index, -delta);
            refinementData.crs_fc.computeDensity(refinementData.fc);
            llk2 = sigmaAMinimize.calculateLikelihood();
            fd = (llk1 - llk2) / (2.0 * delta);
            System.out.print(String.format(" Z A: %16.8f FD: %16.8f Error: %16.8f\n", gxyz[2], fd, gxyz[2] - fd));
            refinementData.crs_fc.deltaZ(index, 0.0);

            nmean++;
            mean += (gxyz[2] / fd - mean) / nmean;

            if (atom.getAnisou(null) == null) {
                double b = atom.getTempFactor();
                atom.setTempFactor(b + delta);
                refinementData.crs_fc.computeDensity(refinementData.fc);
                llk1 = sigmaAMinimize.calculateLikelihood();
                atom.setTempFactor(b - delta);
                refinementData.crs_fc.computeDensity(refinementData.fc);
                llk2 = sigmaAMinimize.calculateLikelihood();
                fd = (llk1 - llk2) / (2.0 * delta);
                System.out.print(String.format(" B A: %16.8f FD: %16.8f Error: %16.8f\n",
                        bg, fd, bg - fd));
                atom.setTempFactor(b);
                nmean++;
                mean += (bg / fd - mean) / nmean;
            } else {
                double anisou[] = atom.getAnisou(null);
                double anisouG[] = atom.getAnisouGradient(null);
                for (int j = 0; j < 6; j++) {
                    double tmpu = anisou[j];
                    anisou[j] = tmpu + b2u(delta);
                    atom.setAnisou(anisou);
                    refinementData.crs_fc.computeDensity(refinementData.fc);
                    llk1 = sigmaAMinimize.calculateLikelihood();
                    anisou[j] = tmpu - b2u(delta);
                    atom.setAnisou(anisou);
                    refinementData.crs_fc.computeDensity(refinementData.fc);
                    llk2 = sigmaAMinimize.calculateLikelihood();
                    fd = (llk1 - llk2) / (2.0 * b2u(delta));
                    atom.getAnisouGradient(anisouG);
                    System.out.print(String.format(" %d A: %16.8f FD: %16.8f Error: %16.8f\n",
                            j, anisouG[j], fd, anisouG[j] - fd));
                    anisou[j] = tmpu;
                    atom.setAnisou(anisou);
                    nmean++;
                    mean += (anisouG[j] / fd - mean) / nmean;
                }
            }
        }

        System.out.println(" Mean df/fd ratio (" + nmean + "): " + mean);
        assertEquals(" Gradient: finite-difference ratio should be approx. 1.0", 1.0, mean, 0.01);
    }
}
