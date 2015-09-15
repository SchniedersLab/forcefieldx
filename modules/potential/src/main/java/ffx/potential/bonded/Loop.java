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
import java.util.ArrayList;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.MolecularAssembly;

import static ffx.numerics.VectorMath.dot;

/**
 * @author Mallory R. Tollefson
 */
public class Loop {

    private static final Logger logger = Logger.getLogger(Loop.class.getName());

    int max_soln = 16;
    double[][] r_n = new double[5][3];
    double[][] r_a = new double[5][3];
    double[][] r_c = new double[5][3];
    public final LoopClosure loopClosure;
    public final SturmMethod sturmMethod;

    public Loop(MolecularAssembly molAss, int firstResidue, int endResidue, boolean writeFile) {

        loopClosure = new LoopClosure();
        sturmMethod = new SturmMethod();

        ArrayList<Atom> atoms = molAss.getBackBoneAtoms();
        boolean bool1 = true;
        int i = 0;

        while (bool1) {
            Atom atom = atoms.get(i);
            int resID = atom.getResidueNumber();

            if (resID > endResidue) {
                //terminates the collection of atom coordinates
                bool1 = false;
            }

            if (resID < firstResidue) {
                i++;
            } else if (resID >= firstResidue && resID <= endResidue) {
                int ir = resID - firstResidue;
                String atmname = atom.getAtomType().name;
                String ca = "CA";
                String c = "C";
                String n = "N";
                //coordinateArray temporarily holds the coordinates of a
                //      specific atom from the pdb file
                double[] coordinateArray = atom.getXYZ(null);

                if (atmname.contentEquals(n)) {
                    //Backbone nitrogen coordinates are stored in r_n[]
                    r_n[ir][0] = coordinateArray[0];
                    r_n[ir][1] = coordinateArray[1];
                    r_n[ir][2] = coordinateArray[2];
                } else if (atmname.contentEquals(ca)) {
                    //Backbone alpha carbon coordinates are stored in r_a[]
                    r_a[ir][0] = coordinateArray[0];
                    r_a[ir][1] = coordinateArray[1];
                    r_a[ir][2] = coordinateArray[2];
                } else if (atmname.contentEquals(c)) {
                    //Backbone carbon coordinates are stored in r_c[]
                    r_c[ir][0] = coordinateArray[0];
                    r_c[ir][1] = coordinateArray[1];
                    r_c[ir][2] = coordinateArray[2];
                }
                i++;
            }
        }

        /**
         * Method that solves the 16th degree, 3 peptide polynomial.
         */
        double[][][] r_soln_n = new double[max_soln][3][3];
        double[][][] r_soln_a = new double[max_soln][3][3];
        double[][][] r_soln_c = new double[max_soln][3][3];
        int[] n_soln = new int[1];
         
        loopClosure.solve3PepPoly(r_n[1], r_a[1], r_a[3], r_c[3], r_soln_n, r_soln_a, r_soln_c, n_soln);

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Starting res.:             %d\n", firstResidue));
        sb.append(String.format(" Ending res.:               %d\n", endResidue));
        sb.append(String.format(" No. of solutions:          %d\n", n_soln[0]));
        logger.info(sb.toString());

        int counter = 0;
        double[] dr = new double[3];

        for (int k = 0; k < n_soln[0]; k++) {
            for (i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    r_n[i + 1][j] = r_soln_n[k][i][j];
                    r_a[i + 1][j] = r_soln_a[k][i][j];
                    r_c[i + 1][j] = r_soln_c[k][i][j];
                }
            }
            double sum = 0.0;
           
            for (i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    dr[j] = r_soln_n[k][i][j] - r_n[i + 1][j];
                }
                sum += dot(dr, dr);
                for (int j = 0; j < 3; j++) {
                    dr[j] = r_soln_a[k][i][j] - r_a[i + 1][j];
                }
                sum += dot(dr, dr);
                for (int j = 0; j < 3; j++) {
                    dr[j] = r_soln_c[k][i][j] - r_c[i + 1][j];
                }
                sum += dot(dr, dr);
            }

            //Here the sum is now equal to the rmsd.
            sum = sqrt(sum / 9.0);

            StringBuilder string = new StringBuilder();
            string.append(String.format("Rmsd for solution #" + (k + 1) + " is " + sum + "\n"));
            logger.info(string.toString());
            counter++;
            File fileName = sturmMethod.writePDBBackbone(r_n, r_a, r_c,
                    firstResidue, endResidue, molAss, counter, writeFile);
            string.append(String.format("Recording the solution #" + (k + 1) + " in " + fileName + ".\n"));
        }
    }

    // Only for JUnit testing.
    public double[][] getRN() {
        return r_n;
    }

    public double[][] getRA() {
        return r_a;
    }

    public double[][] getRC() {
        return r_c;
    }
}
