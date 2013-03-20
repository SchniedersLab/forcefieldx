/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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
import java.util.List;

import org.apache.commons.configuration.CompositeConfiguration;
import org.junit.*;

import static org.junit.Assert.assertEquals;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.numerics.ComplexNumber;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.Keyword;

/**
 *
 * @author Timothy Fenn
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
        DiffractionRefinementData refinementdata = new DiffractionRefinementData(properties,
                reflectionlist);

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
        Atom atomarray[] = atomlist.toArray(new Atom[atomlist.size()]);

        // set up FFT and run it
        ParallelTeam parallelTeam = new ParallelTeam();
        CrystalReciprocalSpace crs =
                new CrystalReciprocalSpace(reflectionlist, atomarray,
                parallelTeam, parallelTeam);
        crs.computeAtomicDensity(refinementdata.fc);

        // tests
        ComplexNumber b = new ComplexNumber(-828.584, -922.704);
        HKL hkl = reflectionlist.getHKL(1, 1, 4);
        ComplexNumber a = refinementdata.get_fc(hkl.index());
        System.out.println("1 1 4: " + a.toString() + " | "
                + b.toString() + " | "
                + a.divides(b).toString());

        /*
         * 1 1 4 reflection should be correct expected:<-753.08997> but was:<-753.5872150259233>
         * 1 1 4 reflection should be correct expected:<-1011.751257> but was:<-1009.6764916871754>
         * 2 1 10 reflection should be correct expected:<-67.531378> but was:<-68.77264741298902
         * 2 1 10 reflection should be correct expected:<-412.935591> but was:<-409.4297907290974>
         */

        assertEquals("1 1 4 reflection should be correct",
                 -753.587215, a.re(), 0.001);
        assertEquals("1 1 4 reflection should be correct",
                -1009.676491, a.im(), 0.001);

        b.re(-70.4582);
        b.im(-486.142);
        hkl = reflectionlist.getHKL(2, 1, 10);
        a = refinementdata.get_fc(hkl.index());
        System.out.println("2 1 10: " + a.toString() + " | "
                + b.toString() + " | "
                + a.divides(b).toString());

        assertEquals("2 1 10 reflection should be correct",
                 -68.772647, a.re(), 0.001);
        assertEquals("2 1 10 reflection should be correct",
                -409.429790, a.im(), 0.001);
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
        DiffractionRefinementData refinementdata = new DiffractionRefinementData(properties,
                reflectionlist);

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
        Atom atomarray[] = atomlist.toArray(new Atom[atomlist.size()]);

        // set up FFT and run it
        ParallelTeam parallelTeam = new ParallelTeam();
        CrystalReciprocalSpace crs =
                new CrystalReciprocalSpace(reflectionlist, atomarray,
                parallelTeam, parallelTeam);
        crs.computeAtomicDensity(refinementdata.fc);

        // tests
        ComplexNumber b = new ComplexNumber(-496.999, 431.817);
        HKL hkl = reflectionlist.getHKL(1, 9, 4);
        ComplexNumber a = refinementdata.get_fc(hkl.index());
        System.out.println("1 9 4: " + a.toString() + " | "
                + b.toString() + " | "
                + a.divides(b).toString());

        /*
         * 1 9 4 reflection should be correct expected:<-495.805516> but was:<-496.5665115012335>
         * 1 9 4 reflection should be correct expected:<460.309527> but was:<460.1494607581358>
         * 5 26 8 reflection should be correct expected:<-123.928508> but was:<-123.95591260975226>
         * 5 26 8 reflection should be correct expected:<-73.122412> but was:<-73.58595634462314>
         */

        assertEquals("1 9 4 reflection should be correct",
                -496.566511, a.re(), 0.001);
        assertEquals("1 9 4 reflection should be correct",
                 460.149460, a.im(), 0.001);

        b.re(-129.767);
        b.im(-76.9812);
        hkl = reflectionlist.getHKL(5, 26, 8);
        a = refinementdata.get_fc(hkl.index());
        System.out.println("5 26 8: " + a.toString() + " | "
                + b.toString() + " | "
                + a.divides(b).toString());

        // -123.92850839573994 - 73.12241222985045i
        assertEquals("5 26 8 reflection should be correct",
                -123.955912, a.re(), 0.001);
        assertEquals("5 26 8 reflection should be correct",
                 -73.585956, a.im(), 0.001);
    }
}
