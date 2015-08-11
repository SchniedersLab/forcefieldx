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
import java.util.logging.Logger;

import ffx.numerics.VectorMath;
import ffx.potential.MolecularAssembly;

/**
 * @author Mallory R. Tollefson
 */
public class Loop {

    private Logger logger = Logger.getLogger(Loop.class.getName());

    private static final LoopClosure loopClosure = new LoopClosure();
    public static final SturmMethod sturmMethod = new SturmMethod();

    int max_soln = 16;
    char[] out_pdb = new char[100];
    int i, j, k;
    int[] n_soln = new int[1];
    double[][] r_n = new double[5][3];
    double[][] r_a = new double[5][3];  
    double[][] r_c = new double[5][3];
    double[][][] r_soln_n = new double[max_soln][3][3];
    double[][][] r_soln_a = new double[max_soln][3][3];
    double[][][] r_soln_c = new double[max_soln][3][3];
    double sum;
    double[] dr = new double[3];
    MolecularAssembly molAss;
    boolean bool = false;
    int counter = 0;

    public Loop(MolecularAssembly molAss, int stt_res, int end_res, boolean writeFile) {
        this.molAss = molAss;

        loopClosure.initializeLoopClosure();

        for (i = 0; i < 1; i++) {
            int res_no;
            int ir;
            boolean bool1 = true;
            i = 0;

            while (bool1) {
                res_no = molAss.getBackBoneAtoms().get(i).getResidueNumber();
                String name = molAss.getBackBoneAtoms().get(i).getName();
                String resName = molAss.getBackBoneAtoms().get(i).getResidueName();

                if (res_no > end_res) {
                    //terminates the collection of atom coordinates
                    bool1 = false;
                }

                if (res_no < stt_res) {
                    i++;
                } else if (res_no >= stt_res && res_no <= end_res) {
                    ir = res_no - stt_res;

                    String atmname;
                    atmname = molAss.getBackBoneAtoms().get(i).getAtomType().name;

                    String ca = "CA";
                    String c = "C";
                    String n = "N";

                    //coordinateArray temporarily holds the coordinates of a
                    //      specific atom from the pdb file
                    double[] coordinateArray = molAss.getBackBoneAtoms().get(i).getXYZ();

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

            //Method that solves the 16th degree, 3 peptide polynomial
            loopClosure.solve3PepPoly(r_n[1], r_a[1], r_a[3], r_c[3], r_soln_n, r_soln_a, r_soln_c, n_soln);

            StringBuilder sb = new StringBuilder();
            sb.append(String.format("Starting res.:             %d\n", stt_res));
            sb.append(String.format("Ending res.:               %d\n", end_res));
            sb.append(String.format("No. of solutions:          %d\n", n_soln[0]));
            logger.info(sb.toString());

                for (k = 0; k < n_soln[0]; k++) {
                    for (i = 0; i < 3; i++) {
                        for (j = 0; j < 3; j++) {
                            r_n[i + 1][j] = r_soln_n[k][i][j];
                            r_a[i + 1][j] = r_soln_a[k][i][j];
                            r_c[i + 1][j] = r_soln_c[k][i][j];
                        }
                    }
                    sum = 0.0e0;

                    for (i = 0; i < 3; i++) {
                        for (j = 0; j < 3; j++) {
                            dr[j] = r_soln_n[k][i][j] - r_n[i+1][j];
                        }
                        sum += VectorMath.dot(dr, dr);
                        for (j = 0; j < 3; j++) {
                            dr[j] = r_soln_a[k][i][j] - r_a[i+1][j];
                        }
                        sum += VectorMath.dot(dr, dr);
                        for (j = 0; j < 3; j++) {
                            dr[j] = r_soln_c[k][i][j] - r_c[i+1][j];
                        }
                        sum += VectorMath.dot(dr, dr);
                    }
                    
                    //Here the sum is now equal to the rmsd.
                    sum = Math.sqrt(sum / 9.0e0);

                    StringBuilder string = new StringBuilder();
                    string.append(String.format("Rmsd for solution #" + (k + 1) + " is " + sum + "\n"));
                    logger.info(string.toString());

                    counter++;

                    File fileName = sturmMethod.writePDBBackbone(r_n, r_a, r_c, stt_res, end_res, molAss, counter, writeFile);
                    StringBuilder string1 = new StringBuilder();
                    string.append(String.format("Recording the solution #" + (k + 1) + " in " + fileName + ".\n"));
                }
            
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
