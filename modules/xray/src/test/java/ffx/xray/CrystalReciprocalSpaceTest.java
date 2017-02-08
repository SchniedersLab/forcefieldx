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
import java.util.List;

import org.apache.commons.configuration.CompositeConfiguration;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.numerics.ComplexNumber;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.utils.PotentialsUtils;
import ffx.utilities.Keyword;

/**
 * @author Timothy D. Fenn
 */
public class CrystalReciprocalSpaceTest {

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
        PotentialsUtils potutil = new PotentialsUtils();
        MolecularAssembly mola = potutil.open(structure);
        CompositeConfiguration properties = mola.getProperties();

        Crystal crystal = new Crystal(39.767, 51.750, 132.938,
                90.00, 90.00, 90.00, "P212121");
        Resolution resolution = new Resolution(1.45);

        ReflectionList reflectionList = new ReflectionList(crystal, resolution);
        DiffractionRefinementData refinementData = new DiffractionRefinementData(properties,
                reflectionList);

        mola.finalize(true, mola.getForceField());
        ForceFieldEnergy energy = mola.getPotentialEnergy();

        List<Atom> atomList = mola.getAtomList();
        Atom atomArray[] = atomList.toArray(new Atom[atomList.size()]);

        // set up FFT and run it
        ParallelTeam parallelTeam = new ParallelTeam();
        CrystalReciprocalSpace crs
                = new CrystalReciprocalSpace(reflectionList, atomArray,
                        parallelTeam, parallelTeam);
        crs.computeAtomicDensity(refinementData.fc);

        // tests
        ComplexNumber b = new ComplexNumber(-828.584, -922.704);
        HKL hkl = reflectionList.getHKL(1, 1, 4);
        ComplexNumber a = refinementData.getFc(hkl.index());
        System.out.println("1 1 4: " + a.toString() + " | "
                + b.toString() + " | "
                + a.divides(b).toString());

        assertEquals("1 1 4 reflection should be correct",
                -753.4722104328416, a.re(), 0.0001);
        assertEquals("1 1 4 reflection should be correct",
                -1012.1341308707799, a.im(), 0.0001);

        b.re(-70.4582);
        b.im(-486.142);
        hkl = reflectionList.getHKL(2, 1, 10);
        a = refinementData.getFc(hkl.index());
        System.out.println("2 1 10: " + a.toString() + " | "
                + b.toString() + " | "
                + a.divides(b).toString());

        assertEquals("2 1 10 reflection should be correct",
                -69.39660884054359, a.re(), 0.0001);
        assertEquals("2 1 10 reflection should be correct",
                -412.0147625765328, a.im(), 0.0001);
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

        ReflectionList reflectionList = new ReflectionList(crystal, resolution);
        DiffractionRefinementData refinementData = new DiffractionRefinementData(properties,
                reflectionList);

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
        Atom atomArray[] = atomList.toArray(new Atom[atomList.size()]);

        // set up FFT and run it
        ParallelTeam parallelTeam = new ParallelTeam();
        CrystalReciprocalSpace crs
                = new CrystalReciprocalSpace(reflectionList, atomArray,
                        parallelTeam, parallelTeam);
        crs.computeAtomicDensity(refinementData.fc);

        // tests
        ComplexNumber b = new ComplexNumber(-496.999, 431.817);
        HKL hkl = reflectionList.getHKL(1, 9, 4);
        ComplexNumber a = refinementData.getFc(hkl.index());
        System.out.println("1 9 4: " + a.toString() + " | "
                + b.toString() + " | "
                + a.divides(b).toString());

        assertEquals("1 9 4 reflection should be correct",
                -493.7799429881329, a.re(), 0.0001);
        assertEquals("1 9 4 reflection should be correct",
                460.7022632345927, a.im(), 0.0001);

        b.re(-129.767);
        b.im(-76.9812);
        hkl = reflectionList.getHKL(5, 26, 8);
        a = refinementData.getFc(hkl.index());
        System.out.println("5 26 8: " + a.toString() + " | "
                + b.toString() + " | "
                + a.divides(b).toString());

        assertEquals("5 26 8 reflection should be correct",
                -123.05535567943377, a.re(), 0.0001);
        assertEquals("5 26 8 reflection should be correct",
                -74.59007322382718, a.im(), 0.0001);
    }
}
