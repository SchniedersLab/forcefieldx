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

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import edu.rit.pj.ParallelTeam;
import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.Keyword;

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
public class ScaleBulkMinimizeTest {

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                    {true,
                        "NSF D2 domain test",
                        "ffx/xray/structures/1NSF.pdb",
                        "ffx/xray/structures/1NSF.mtz",
                        null,
                        25.39,
                        25.54,
                        0.8940,
                        0.1515},
                    {true,
                        "SNARE complex",
                        "ffx/xray/structures/1N7S.pdb",
                        "ffx/xray/structures/1N7S.mtz",
                        null,
                        19.60,
                        21.76,
                        0.9311,
                        0.1361}
                });
    }
    private final String info;
    private final String pdbname;
    private final String mtzname;
    private final String cifname;
    private final CrystalStats crystalstats;
    private final double r;
    private final double rfree;
    private final double sigmaa;
    private final double sigmaw;
    private final boolean ci;
    private final boolean ciOnly;

    public ScaleBulkMinimizeTest(boolean ciOnly,
            String info, String pdbname, String mtzname, String cifname,
            double r, double rfree, double sigmaa, double sigmaw) {
        this.ciOnly = ciOnly;
        this.info = info;
        this.pdbname = pdbname;
        this.mtzname = mtzname;
        this.cifname = cifname;
        this.r = r;
        this.rfree = rfree;
        this.sigmaa = sigmaa;
        this.sigmaw = sigmaw;

        String ffxCi = System.getProperty("ffx.ci");
        if (ffxCi != null && ffxCi.equalsIgnoreCase("true")) {
            ci = true;
        } else {
            ci = false;
        }

        if (!ci && ciOnly) {
            crystalstats = null;
            return;
        }

        int index = pdbname.lastIndexOf(".");
        String name = pdbname.substring(0, index);

        // load the structure
        ClassLoader cl = this.getClass().getClassLoader();
        File structure = new File(cl.getResource(pdbname).getPath());
        File mtzfile = null, ciffile = null;
        if (mtzname != null) {
            mtzfile = new File(cl.getResource(mtzname).getPath());
        } else {
            ciffile = new File(cl.getResource(cifname).getPath());
        }

        // load any properties associated with it
        CompositeConfiguration properties = Keyword.loadProperties(structure);

        // read in Fo/sigFo/FreeR
        MTZFilter mtzfilter = new MTZFilter();
        CIFFilter ciffilter = new CIFFilter();
        ReflectionList reflectionlist;
        Crystal crystal = Crystal.checkProperties(properties);
        Resolution resolution = Resolution.checkProperties(properties);
        if (crystal == null || resolution == null) {
            if (mtzname != null) {
                reflectionlist = mtzfilter.getReflectionList(mtzfile);
            } else {
                reflectionlist = ciffilter.getReflectionList(ciffile);
            }
        } else {
            reflectionlist = new ReflectionList(crystal, resolution);
        }

        RefinementData refinementdata = new RefinementData(properties,
                reflectionlist);
        if (mtzname != null) {
            assertTrue(info + " mtz file should be read in without errors",
                    mtzfilter.readFile(mtzfile, reflectionlist, refinementdata));
        } else {
            assertTrue(info + " cif file should be read in without errors",
                    ciffilter.readFile(ciffile, reflectionlist, refinementdata));
        }

        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
        ForceField forceField = forceFieldFilter.parse();

        // associate molecular assembly with the structure, set up forcefield
        MolecularAssembly molecularAssembly = new MolecularAssembly(name);
        molecularAssembly.setFile(structure);
        molecularAssembly.setForceField(forceField);
        PDBFilter pdbfile = new PDBFilter(molecularAssembly);
        pdbfile.setForceField(forceField);
        pdbfile.setProperties(properties);
        pdbfile.readFile();

        List<Atom> atomlist = molecularAssembly.getAtomList();
        Atom atomarray[] = atomlist.toArray(new Atom[atomlist.size()]);

        // set up FFT and run it
        ParallelTeam parallelTeam = new ParallelTeam();
        CrystalReciprocalSpace crs = new CrystalReciprocalSpace(reflectionlist,
                atomarray, parallelTeam, parallelTeam, false);
        crs.computeDensity(refinementdata.fc);
        refinementdata.setCrystalReciprocalSpaceFc(crs);
        crs = new CrystalReciprocalSpace(reflectionlist, atomarray,
                parallelTeam, parallelTeam, true);
        crs.computeDensity(refinementdata.fs);
        refinementdata.setCrystalReciprocalSpaceFs(crs);

        /*
        try {
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("/tmp/foo.cns")));
        out.println("ANOMalous=FALSE");
        out.println("DECLare NAME=FC DOMAin=RECIprocal TYPE=COMP END");
        for (HKL ih : reflectionlist.hkllist) {
        if (ih.h() == 0
        && ih.k() == 0
        && ih.l() == 0) {
        continue;
        }
        double fc[] = refinementdata.fs[ih.index()];
        out.printf("INDE %5d%5d%5d FC= ", ih.h(), ih.k(), ih.l());
        out.printf("%10.3f%10.3f\n",
        Math.hypot(fc[0], fc[1]),
        Math.toDegrees(Math.atan2(fc[1], fc[0])));
        }
        out.close();
        } catch (Exception e) {
        System.out.println("error: " + e.getMessage());
        }
         */

        ScaleBulkMinimize scalebulkminimize =
                new ScaleBulkMinimize(reflectionlist, refinementdata, crs);
        /*
        if (refinementdata.solvent_n > 1) {
            scalebulkminimize.minimize(7, 1e-2);
            scalebulkminimize.GridOptimize();
        }
        */
        scalebulkminimize.minimize(7, 1e-4);

        SigmaAMinimize sigmaaminimize = new SigmaAMinimize(reflectionlist,
                refinementdata);
        sigmaaminimize.minimize(7, 1e-1);

        SplineMinimize splineminimize = new SplineMinimize(reflectionlist,
                refinementdata, refinementdata.spline,
                SplineEnergy.Type.FOFC);
        splineminimize.minimize(7, 1e-5);

        crystalstats = new CrystalStats(reflectionlist, refinementdata);
    }

    @Test
    public void testCrystalStats() {
        if (!ci && ciOnly) {
            return;
        }

        crystalstats.print_scalestats();
        crystalstats.print_hklstats();
        crystalstats.print_snstats();
        crystalstats.print_rstats();

        assertEquals(info + " R value should be correct",
                r, crystalstats.get_r(), 0.01);
        assertEquals(info + " Rfree value should be correct",
                rfree, crystalstats.get_rfree(), 0.01);
        assertEquals(info + " sigmaA s value should be correct",
                sigmaa, crystalstats.get_sigmaa(), 0.01);
        assertEquals(info + " sigmaA w value should be correct",
                sigmaw, crystalstats.get_sigmaw(), 0.01);
    }
}
