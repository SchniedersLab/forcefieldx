/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.xray;

import static ffx.numerics.VectorMath.b2u;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import edu.rit.pj.ParallelTeam;
import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.Keyword;
import ffx.xray.CrystalReciprocalSpace.SolventModel;
import ffx.xray.RefinementMinimize.RefinementMode;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import static org.junit.Assert.*;

/**
 *
 * @author fennt
 */
@RunWith(Parameterized.class)
public class FiniteDifferenceTest {

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                    {false,
                        "ala met anisou",
                        SolventModel.NONE,
                        new int[]{91, 105, 119},
                        "ffx/xray/structures/alamet.pdb",
                        "ffx/xray/structures/alamet.mtz"}
                });
    }
    private final String info;
    private final int solventmodel;
    private final String pdbname;
    private final String mtzname;
    private final boolean ci;
    private final boolean ciOnly;
    private final Atom atomarray[];
    private final int atoms[];
    private final DiffractionRefinementData refinementdata;
    private final SigmaAMinimize sigmaaminimize;

    public FiniteDifferenceTest(boolean ciOnly,
            String info, int solventmodel, int[] atoms,
            String pdbname, String mtzname) {
        this.ciOnly = ciOnly;
        this.info = info;
        this.solventmodel = solventmodel;
        this.atoms = atoms;
        this.pdbname = pdbname;
        this.mtzname = mtzname;

        String ffxCi = System.getProperty("ffx.ci");
        if (ffxCi != null && ffxCi.equalsIgnoreCase("true")) {
            ci = true;
        } else {
            ci = false;
        }

        if (!ci && ciOnly) {
            atomarray = null;
            refinementdata = null;
            sigmaaminimize = null;
            return;
        }

        int index = pdbname.lastIndexOf(".");
        String name = pdbname.substring(0, index);

        // load the structure
        ClassLoader cl = this.getClass().getClassLoader();
        File structure = new File(cl.getResource(pdbname).getPath());
        File mtzfile = null, ciffile = null;
        mtzfile = new File(cl.getResource(mtzname).getPath());

        // load any properties associated with it
        CompositeConfiguration properties = Keyword.loadProperties(structure);

        // read in Fo/sigFo/FreeR
        MTZFilter mtzfilter = new MTZFilter();
        CIFFilter ciffilter = new CIFFilter();
        ReflectionList reflectionlist;
        Crystal crystal = Crystal.checkProperties(properties);
        Resolution resolution = Resolution.checkProperties(properties);
        if (crystal == null || resolution == null) {
            reflectionlist = mtzfilter.getReflectionList(mtzfile);
        } else {
            reflectionlist = new ReflectionList(crystal, resolution);
        }

        refinementdata = new DiffractionRefinementData(properties,
                reflectionlist);
        assertTrue(info + " mtz file should be read in without errors",
                mtzfilter.readFile(mtzfile, reflectionlist, refinementdata,
                properties));

        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
        ForceField forceField = forceFieldFilter.parse();

        // associate molecular assembly with the structure, set up forcefield
        MolecularAssembly molecularAssembly = new MolecularAssembly(name);
        molecularAssembly.setFile(structure);
        molecularAssembly.setForceField(forceField);
        PDBFilter pdbfile = new PDBFilter(structure, molecularAssembly, forceField, properties);
        pdbfile.readFile();
        molecularAssembly.finalize(true);
        ForceFieldEnergy energy = new ForceFieldEnergy(molecularAssembly);

        List<Atom> atomlist = molecularAssembly.getAtomList();
        atomarray = atomlist.toArray(new Atom[atomlist.size()]);
        boolean use_3g = properties.getBoolean("use_3g", true);

        // initialize atomic form factors
        for (int i = 0; i < atomarray.length; i++) {
            XRayFormFactor atomff =
                    new XRayFormFactor(atomarray[i], use_3g, 2.0);
            atomarray[i].setFormFactorIndex(atomff.ffindex);

            if (atomarray[i].getOccupancy() == 0.0) {
                atomarray[i].setFormFactorWidth(1.0);
                continue;
            }

            double arad = 2.4;
            double xyz[] = new double[3];
            xyz[0] = atomarray[i].getX() + arad;
            xyz[1] = atomarray[i].getY();
            xyz[2] = atomarray[i].getZ();
            while (true) {
                double rho = atomff.rho(0.0, 1.0, xyz);
                if (rho > 0.1) {
                    arad += 0.5;
                } else if (rho > 0.001) {
                    arad += 0.1;
                } else {
                    arad += 0.2;
                    atomarray[i].setFormFactorWidth(arad);
                    break;
                }
                xyz[0] = atomarray[i].getX() + arad;
            }
        }

        // set up FFT and run it
        ParallelTeam parallelTeam = new ParallelTeam();
        CrystalReciprocalSpace crs = new CrystalReciprocalSpace(reflectionlist,
                atomarray, parallelTeam, parallelTeam, false, false);
        crs.computeDensity(refinementdata.fc);
        refinementdata.setCrystalReciprocalSpace_fc(crs);
        crs = new CrystalReciprocalSpace(reflectionlist,
                atomarray, parallelTeam, parallelTeam, true, false, solventmodel);
        crs.computeDensity(refinementdata.fs);
        refinementdata.setCrystalReciprocalSpace_fs(crs);

        ScaleBulkMinimize scalebulkminimize =
                new ScaleBulkMinimize(reflectionlist, refinementdata, crs);
        scalebulkminimize.minimize(7, 1e-4);

        sigmaaminimize = new SigmaAMinimize(reflectionlist,
                refinementdata);
        sigmaaminimize.minimize(7, 1.0);

        SplineMinimize splineminimize = new SplineMinimize(reflectionlist,
                refinementdata, refinementdata.spline,
                SplineEnergy.Type.FOFC);
        splineminimize.minimize(7, 1e-5);

        CrystalStats crystalstats = new CrystalStats(reflectionlist, refinementdata);
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
        refinementdata.crs_fc.computeAtomicGradients(refinementdata.dfc,
                refinementdata.freer, refinementdata.rfreeflag,
                RefinementMode.COORDINATES_AND_BFACTORS);

        int natoms = atomarray.length;
        double llk0 = refinementdata.llkr;
        double llk1, llk2, fd;
        double mean = 0.0, nmean = 0.0;
        double gxyz[] = new double[3];
        for (int i = 0; i < atoms.length; i++) {
            Atom atom = atomarray[atoms[i]];
            int index = atom.getXYZIndex() - 1;
            if (atom.getOccupancy() == 0.0) {
                continue;
            }

            System.out.println("atom " + i + ": (" + atom.getXYZIndex() + ") " + atom.toString());
            atom.getXYZGradient(gxyz);
            double bg = atom.getTempFactorGradient();
            double anisoug[] = null;
            if (atom.getAnisou() != null) {
                anisoug = atom.getAnisouGradient();
            }

            refinementdata.crs_fc.deltaX(index, delta);
            refinementdata.crs_fc.computeDensity(refinementdata.fc);
            llk1 = sigmaaminimize.calculateLikelihood();
            refinementdata.crs_fc.deltaX(index, -delta);
            refinementdata.crs_fc.computeDensity(refinementdata.fc);
            llk2 = sigmaaminimize.calculateLikelihood();
            fd = (llk1 - llk2) / (2.0 * delta);
            System.out.print(String.format("+x: %g -x: %g dfx: %g fdx: %g ratio: %g\n",
                    llk1 - llk0, llk2 - llk0, gxyz[0], fd, gxyz[0] / fd));
            refinementdata.crs_fc.deltaX(index, 0.0);

            nmean++;
            mean += (gxyz[0] / fd - mean) / nmean;

            refinementdata.crs_fc.deltaY(index, delta);
            refinementdata.crs_fc.computeDensity(refinementdata.fc);
            llk1 = sigmaaminimize.calculateLikelihood();
            refinementdata.crs_fc.deltaY(index, -delta);
            refinementdata.crs_fc.computeDensity(refinementdata.fc);
            llk2 = sigmaaminimize.calculateLikelihood();
            fd = (llk1 - llk2) / (2.0 * delta);
            System.out.print(String.format("+y: %g -y: %g dfy: %g fdy: %g ratio: %g\n",
                    llk1 - llk0, llk2 - llk0, gxyz[1], fd, gxyz[1] / fd));
            refinementdata.crs_fc.deltaY(index, 0.0);

            nmean++;
            mean += (gxyz[1] / fd - mean) / nmean;

            refinementdata.crs_fc.deltaZ(index, delta);
            refinementdata.crs_fc.computeDensity(refinementdata.fc);
            llk1 = sigmaaminimize.calculateLikelihood();
            refinementdata.crs_fc.deltaZ(index, -delta);
            refinementdata.crs_fc.computeDensity(refinementdata.fc);
            llk2 = sigmaaminimize.calculateLikelihood();
            fd = (llk1 - llk2) / (2.0 * delta);
            System.out.print(String.format("+z: %g -z: %g dfz: %g fdz: %g ratio: %g\n",
                    llk1 - llk0, llk2 - llk0, gxyz[2], fd, gxyz[2] / fd));
            refinementdata.crs_fc.deltaZ(index, 0.0);

            nmean++;
            mean += (gxyz[2] / fd - mean) / nmean;

            if (atom.getAnisou() == null) {
                double b = atom.getTempFactor();
                atom.setTempFactor(b + delta);
                refinementdata.crs_fc.computeDensity(refinementdata.fc);
                llk1 = sigmaaminimize.calculateLikelihood();
                atom.setTempFactor(b - delta);
                refinementdata.crs_fc.computeDensity(refinementdata.fc);
                llk2 = sigmaaminimize.calculateLikelihood();
                fd = (llk1 - llk2) / (2.0 * delta);
                System.out.print(String.format("+B: %g -B: %g dfB: %g fdB: %g ratio: %g\n",
                        llk1 - llk0, llk2 - llk0, bg, fd, bg / fd));
                atom.setTempFactor(b);

                nmean++;
                mean += (bg / fd - mean) / nmean;
            } else {
                double anisou[] = atom.getAnisou();
                for (int j = 0; j < 6; j++) {
                    double tmpu = anisou[j];
                    anisou[j] = tmpu + b2u(delta);
                    refinementdata.crs_fc.computeDensity(refinementdata.fc);
                    llk1 = sigmaaminimize.calculateLikelihood();
                    anisou[j] = tmpu - b2u(delta);
                    refinementdata.crs_fc.computeDensity(refinementdata.fc);
                    llk2 = sigmaaminimize.calculateLikelihood();
                    fd = (llk1 - llk2) / (2.0 * b2u(delta));
                    System.out.print(String.format("+u%d: %g -u%d: %g dfu%d: %g fdu%d: %g ratio: %g\n",
                            j, llk1 - llk0, j, llk2 - llk0, j, anisoug[j], j, fd, anisoug[j] / fd));
                    anisou[j] = tmpu;

                    nmean++;
                    mean += (anisoug[j] / fd - mean) / nmean;
                }
            }
        }

        System.out.println("mean df/fd ratio (" + nmean + "): " + mean);

        assertEquals("gradient : finite difference ratio should be approx. 1.0",
                1.00, mean, 0.01);
    }
}
