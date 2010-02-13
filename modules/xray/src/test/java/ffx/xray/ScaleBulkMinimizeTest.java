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
import java.util.List;

import edu.rit.pj.ParallelTeam;
import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.numerics.ComplexNumber;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.Keyword;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fennt
 */
public class ScaleBulkMinimizeTest {

    public ScaleBulkMinimizeTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of minimize method, of class ScaleMinimize.
     */
    @Test
    public void testMinimize() {
        String mtzfilename = "ffx/xray/structures/1N7S.mtz";
        String pdbfilename = "ffx/xray/structures/1N7S.pdb";
        // String mtzfilename = "ffx/xray/structures/2DRM.mtz";
        // String pdbfilename = "ffx/xray/structures/2DRM.pdb";
        int index = pdbfilename.lastIndexOf(".");
        String name = pdbfilename.substring(0, index);

        // load the structure
        ClassLoader cl = this.getClass().getClassLoader();
        File structure = new File(cl.getResource(pdbfilename).getPath());
        File mtzfile = new File(cl.getResource(mtzfilename).getPath());

        // read in Fo/sigFo/FreeR
        MTZFilter mtzfilter = new MTZFilter();
        ReflectionList reflectionlist = mtzfilter.getReflectionList(mtzfile);
        RefinementData refinementdata =
                new RefinementData(reflectionlist.hkllist.size());
        assertTrue("mtz file should be read in without errors",
                mtzfilter.readFile(mtzfile, reflectionlist, refinementdata));

        // load any properties associated with it
        CompositeConfiguration properties = Keyword.loadProperties(structure);
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties, null);
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
        ParallelTeam parallelTeam = new ParallelTeam(1);
        CrystalReciprocalSpace crs =
                new CrystalReciprocalSpace(atomarray, atomarray.length,
                parallelTeam, reflectionlist, false);
        crs.permanent(refinementdata.fc);
        crs = new CrystalReciprocalSpace(atomarray, atomarray.length,
                parallelTeam, reflectionlist, true);
        crs.permanent(refinementdata.fs);

        /*
        try {
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("/tmp/foo.cns")));
        out.println("ANOMalous=FALSE");
        out.println("DECLare NAME=FC DOMAin=RECIprocal TYPE=COMP END");
        for (HKL ih : reflectionlist.hkllist) {
        if (ih.allowed() == 0.0) {
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
                new ScaleBulkMinimize(reflectionlist, refinementdata);
        scalebulkminimize.minimize(7, 1e-4);

        SigmaAMinimize sigmaaminimize = new SigmaAMinimize(reflectionlist,
                refinementdata);
        sigmaaminimize.minimize(7, 1.0);

        System.out.println("final sigmaA params:");
        for (int i = 0; i < refinementdata.nparams; i++) {
            System.out.printf("%8g  %8g\n", refinementdata.sigmaa[i],
                    refinementdata.sigmaw[i]);
        }
    }
}
