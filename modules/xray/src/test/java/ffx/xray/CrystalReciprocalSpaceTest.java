/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.xray;

import java.io.File;
import java.util.List;
// temporary
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;

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

import edu.rit.pj.ParallelTeam;

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
    public void test2DRMPermanent() {
        // String filename = "ffx/xray/structures/1N7S.pdb";
        // String filename = "ffx/xray/structures/1EXR.pdb";
        String filename = "ffx/xray/structures/2DRM.pdb";
        int index = filename.lastIndexOf(".");
        String name = filename.substring(0, index);

        /*
        Crystal crystal = new Crystal(39.767, 51.750, 132.938,
        90.00, 90.00, 90.00, "P212121");
        Resolution resolution = new Resolution(1.45);
         */
        /*
        Crystal crystal = new Crystal(25.015, 29.415, 52.761, 89.54, 86.10, 82.39, "P1");
        Resolution resolution = new Resolution(1.0);
         */

        // set up the crystal data
        Crystal crystal =
                new Crystal(29.97, 37.86, 44.51, 90.28, 90.11, 90.64, "P1");
        Resolution resolution = new Resolution(1.35);
        ReflectionList reflectionlist = new ReflectionList(crystal, resolution);
        RefinementData refinementdata =
                new RefinementData(reflectionlist.hkllist.size());

        // load the structure
        ClassLoader cl = this.getClass().getClassLoader();
        File structure = new File(cl.getResource(filename).getPath());

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
                parallelTeam, false, reflectionlist, refinementdata.fc);
        crs.permanent();

        // write out reflections (temporary)
        /*
        try {
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("/tmp/foo.cns")));
            out.println("ANOMalous=FALSE");
            out.println("DECLare NAME=FC DOMAin=RECIprocal TYPE=COMP END");
            for (HKL ih : reflectionlist.hkllist) {
                if (ih.allowed == 0) {
                    continue;
                }
                double fc[] = hkldata[ih.index];
                out.printf("INDE %5d%5d%5d FC= ", ih.h, ih.k, ih.l);
                out.printf("%10.3f%10.3f\n",
                        Math.hypot(fc[0], fc[1]),
                        Math.toDegrees(Math.atan2(fc[1], fc[0])));
            }
            out.close();
        } catch (Exception e) {
            String message = "Fatal exception writing structure factors.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
        }
         */

        // tests
        ComplexNumber b = new ComplexNumber(-14.6045, -45.6596);
        HKL hkl = reflectionlist.getHKL(1, 0, 0);
        ComplexNumber a = new ComplexNumber(refinementdata.fc[hkl.index][0],
                refinementdata.fc[hkl.index][1]);
        System.out.println("1 0 0: " + a.toString() + " "
                + a.abs() + " " + a.phase() + " "
                + a.divides(b).toString());

        b.re(113.862);
        b.im(178.684);
        hkl = reflectionlist.getHKL(2, 0, 0);
        a = new ComplexNumber(refinementdata.fc[hkl.index][0],
                refinementdata.fc[hkl.index][1]);
        System.out.println("2 0 0: " + a.toString() + " "
                + a.abs() + " " + a.phase() + " "
                + a.divides(b).toString());
    }
}
