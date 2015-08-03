/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.potential.bonded;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.logging.Logger;

import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertArrayEquals;

import ffx.potential.MolecularAssembly;
import ffx.potential.utils.PotentialsUtils;

/**
 * @author Mallory R. Tollefson
 */
@RunWith(Parameterized.class)
public class LoopClosureTest {

    private static final Logger logger = Logger.getLogger(LoopClosureTest.class.getName());
    private MolecularAssembly molecularAssembly;
    private File structure;
    private Loop loop;

    double[][] xyz_n_test = new double[3][3];
    double[][] xyz_c_test = new double[3][3];
    double[][] xyz_a_test = new double[3][3];
    double[][] xyz_o_test = new double[3][3];

    @Parameters
    public static Collection<Object[]> data() {

        double xyz_n_test[][] = new double[3][3];
        double xyz_a_test[][] = new double[3][3];
        double xyz_c_test[][] = new double[3][3];
        double xyz_o_test[][] = new double[3][3];
        //residue 1
        xyz_n_test[0][0] = 7.773;
        xyz_n_test[0][1] = -9.71;
        xyz_n_test[0][2] = -7.32;
        xyz_a_test[0][0] = 6.331;
        xyz_a_test[0][1] = -9.839;
        xyz_a_test[0][2] = -7.259;
        xyz_c_test[0][0] = 5.886372894231285;
        xyz_c_test[0][1] = -10.55641925089512;
        xyz_c_test[0][2] = -5.994873283542817;
        xyz_o_test[0][0] = 4.7066623518635335;
        xyz_o_test[0][1] = -10.772063009151791;
        xyz_o_test[0][2] = -5.7426213147669065;
        //residue 2
        xyz_n_test[1][0] = 6.267265566616004;
        xyz_n_test[1][1] = -11.821304411459156;
        xyz_n_test[1][2] = -5.840321341761048;
        xyz_a_test[1][0] = 5.873570174757412;
        xyz_a_test[1][1] = -12.55730694668949;
        xyz_a_test[1][2] = -4.654655197113309;
        xyz_c_test[1][0] = 4.522327673119161;
        xyz_c_test[1][1] = -13.229312851952344;
        xyz_c_test[1][2] = -4.836181407502477;
        xyz_o_test[1][0] = 4.000608042124467;
        xyz_o_test[1][1] = -13.903386108027433;
        xyz_o_test[1][2] = -3.955679208709603;
        //residue 3
        xyz_n_test[2][0] = 3.9321584783416297;
        xyz_n_test[2][1] = -13.724556941299172;
        xyz_n_test[2][2] = -3.7520533645561343;
        xyz_a_test[2][0] = 2.642;
        xyz_a_test[2][1] = -14.377;
        xyz_a_test[2][2] = -3.863;
        xyz_c_test[2][0] = 1.658;
        xyz_c_test[2][1] = -13.856;
        xyz_c_test[2][2] = -2.821;
        xyz_o_test[2][0] = 0.5084362754396345;
        xyz_o_test[2][1] = -14.272368583539699;
        xyz_o_test[2][2] = -2.7373896189696216;

        return Arrays.asList(
                new Object[][]{
                    {xyz_n_test, xyz_a_test, xyz_c_test, xyz_o_test}, // constructor arguments for test set 1
                });
    }

    public LoopClosureTest(double[][] xyz_n_test, double[][] xyz_a_test, double[][] xyz_c_test, double[][] xyz_o_test) {
        int stt_res = 1;
        int end_res = 5;
        boolean writeFile = false;
        ClassLoader cl = this.getClass().getClassLoader();
        structure = new File(cl.getResource("ffx/potential/structures/LoopClosureTest.pdb").getPath());
        PotentialsUtils potentialUtils = new PotentialsUtils();
        molecularAssembly = potentialUtils.open(structure.getAbsolutePath())[0];
        loop = new Loop(molecularAssembly, stt_res, end_res, writeFile);

        this.xyz_n_test = xyz_n_test;
        this.xyz_a_test = xyz_a_test;
        this.xyz_c_test = xyz_c_test;
        this.xyz_o_test = xyz_o_test;
    }

    @Before
    public void setUp() {
//        logger.info("\n\n\n\n\n\n\n\n\n\n\nHi, this is Mallory's Test!");
//        int stt_res = 1;
//        int end_res = 5;
//        ClassLoader cl = this.getClass().getClassLoader();
//        structure = new File(cl.getResource("ffx/potential/structures/loopClosureTestPDB.pdb").getPath());
//        PotentialsUtils potentialUtils = new PotentialsUtils();
//        molecularAssembly = potentialUtils.open(structure.getAbsolutePath())[0];
//        Loop loop = new Loop(molecularAssembly, stt_res, end_res);
    }

    @Test
    public void loopTest() {
        double[][] r_a = new double[5][3];
        double[][] r_c = new double[5][3];
        double[][] r_n = new double[5][3];
        double[][] r_o = new double[5][3];

        r_a = loop.getr_a();
        r_c = loop.getr_c();
        r_n = loop.getr_n();
        r_o = Loop.sturmMethod.getr_o();

        //System.out.println("R_O:\n");
        //System.out.println(r_o+"\n\n\n\n");
        //System.out.println(r_o);
        int j = 0;

        for (j = 0; j < 3; j++) {
            assertArrayEquals(r_a[j + 1], xyz_a_test[j], 1e-8);
        }

        for (j = 0; j < 3; j++) {
            assertArrayEquals(r_c[j + 1], xyz_c_test[j], 1e-8);
        }

        for (j = 0; j < 3; j++) {
            assertArrayEquals(r_n[j + 1], xyz_n_test[j], 1e-8);
        }

        for (j = 0; j < 3; j++) {
            assertArrayEquals(r_o[j + 1], xyz_o_test[j], 1e-8);
        }

    }
}
