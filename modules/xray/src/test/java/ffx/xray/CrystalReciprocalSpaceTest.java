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

import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.Keyword;

import edu.rit.pj.ParallelTeam;
import ffx.crystal.HKL;
import ffx.numerics.ComplexNumber;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;

import org.apache.commons.configuration.CompositeConfiguration;
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
public class CrystalReciprocalSpaceTest {

    public CrystalReciprocalSpaceTest() {
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
     * Test of permanent method, of class CrystalReciprocalSpace.
     */
    @Test
    public void test1N7SPermanent() {
        String filename = "ffx/xray/structures/1N7S.pdb";
        int index = filename.lastIndexOf(".");
        String name = filename.substring(0, index);

        // load the structure
        ClassLoader cl = this.getClass().getClassLoader();
        File structure = new File(cl.getResource(filename).getPath());

        // load any properties associated with it
        CompositeConfiguration properties = Keyword.loadProperties(structure);

        Crystal crystal = new Crystal(39.767, 51.750, 132.938,
                90.00, 90.00, 90.00, "P212121");
        Resolution resolution = new Resolution(1.45);

        ReflectionList reflectionlist = new ReflectionList(crystal, resolution);
        RefinementData refinementdata = new RefinementData(properties,
                reflectionlist);

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

        // tests
        ComplexNumber b = new ComplexNumber(-828.584, -922.704);
        HKL hkl = reflectionlist.getHKL(1, 1, 4);
        ComplexNumber a = new ComplexNumber(refinementdata.fc[hkl.index()][0],
                refinementdata.fc[hkl.index()][1]);
        System.out.println("1 1 4: " + a.toString() + " | "
                + b.toString() + " | "
                + a.divides(b).toString());

        assertEquals("1 1 4 reflection should be correct",
                -753.1078686943852, a.re(), 0.001);
        assertEquals("1 1 4 reflection should be correct",
                -991.3980428343319, a.im(), 0.001);

        b.re(-70.4582);
        b.im(-486.142);
        hkl = reflectionlist.getHKL(2, 1, 10);
        a.re(refinementdata.fc[hkl.index()][0]);
        a.im(refinementdata.fc[hkl.index()][1]);
        System.out.println("2 1 10: " + a.toString() + " | "
                + b.toString() + " | "
                + a.divides(b).toString());

        assertEquals("2 1 10 reflection should be correct",
                -80.19700297765023, a.re(), 0.001);
        assertEquals("2 1 10 reflection should be correct",
                -414.60777043854677, a.im(), 0.001);
    }

    @Test
    public void test1NSFPermanent() {
        String filename = "ffx/xray/structures/1NSF.pdb";
        int index = filename.lastIndexOf(".");
        String name = filename.substring(0, index);

        // load the structure
        ClassLoader cl = this.getClass().getClassLoader();
        File structure = new File(cl.getResource(filename).getPath());

        // load any properties associated with it
        CompositeConfiguration properties = Keyword.loadProperties(structure);

        Crystal crystal = new Crystal(115.996, 115.996, 44.13, 90.0, 90.0, 120.0, "P6");
        Resolution resolution = new Resolution(1.89631);

        ReflectionList reflectionlist = new ReflectionList(crystal, resolution);
        RefinementData refinementdata = new RefinementData(properties,
                reflectionlist);

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

        // write out reflections (temporary)
        /*
        try {
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("/tmp/foo.cns")));
        out.println("ANOMalous=FALSE");
        out.println("DECLare NAME=FC DOMAin=RECIprocal TYPE=COMP END");
        for (HKL ih : reflectionlist.hkllist) {
        if (ih.allowed() == 0.0) {
        continue;
        }
        double fc[] = refinementdata.fc[ih.index()];
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

        // tests
        /*
        ComplexNumber b = new ComplexNumber(-828.584, -922.704);
        HKL hkl = reflectionlist.getHKL(1, 1, 4);
        ComplexNumber a = new ComplexNumber(refinementdata.fc[hkl.index()][0],
        refinementdata.fc[hkl.index()][1]);
        System.out.println("1 1 4: " + a.toString() + " | "
        + b.toString() + " | "
        + a.divides(b).toString());

        assertEquals("1 1 4 reflection should be correct",
        -753.1078686943852, a.re(), 0.001);
        assertEquals("1 1 4 reflection should be correct",
        -991.3980428343319, a.im(), 0.001);

        b.re(-70.4582);
        b.im(-486.142);
        hkl = reflectionlist.getHKL(2, 1, 10);
        a.re(refinementdata.fc[hkl.index()][0]);
        a.im(refinementdata.fc[hkl.index()][1]);
        System.out.println("2 1 10: " + a.toString() + " | "
        + b.toString() + " | "
        + a.divides(b).toString());

        assertEquals("2 1 10 reflection should be correct",
        -80.19700297765023, a.re(), 0.001);
        assertEquals("2 1 10 reflection should be correct",
        -414.60777043854677, a.im(), 0.001);
         */
    }
}
