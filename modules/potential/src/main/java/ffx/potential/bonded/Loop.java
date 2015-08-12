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
    char[] in_pdb = new char[100];
    char[] out_prefix = new char[100];
    char[] out_pdb = new char[100];
    String[] res_name = new String[5];
    int n0, i, j, k, n1, n2;
    int[] n_soln = new int[1];
    double[][] r_n = new double[5][3];
    double[][] r_a = new double[5][3];
    double[][] r_c = new double[5][3];
    char[] chain_a = new char[5];
    char[] chain_n = new char[5];
    char[] chain_c = new char[5];
    double[][][] r_soln_n = new double[max_soln][3][3];
    double[][][] r_soln_a = new double[max_soln][3][3];
    double[][][] r_soln_c = new double[max_soln][3][3];
    double rmsd, sum;
    double[] dr = new double[3];
    double[][] r0_n = new double[3][3];
    double[][] r0_a = new double[3][3];
    double[][] r0_c = new double[3][3];
    boolean calc_rmsd;
    int write_out_pdb;
    double[] b_len = new double[6];
    double[] b_ang = new double[7];
    double[] t_ang = new double[2];
    double b_na = 1.45e0, b_ac = 1.52e0, b_cn = 1.33e0;
    double ang_nac = Math.toRadians(111.6), ang_acn = Math.toRadians(117.5);
    double ang_cna = Math.toRadians(120.0);
    MolecularAssembly molAss;
    boolean bool = false;
    int counter = 0;

    public Loop(MolecularAssembly molAss, int stt_res, int end_res, boolean writeFile) {
        this.molAss = molAss;

        //initialize bond lengths based on known constants
        b_len[0] = b_ac;
        b_len[1] = b_cn;
        b_len[2] = b_na;
        b_len[3] = b_ac;
        b_len[4] = b_cn;
        b_len[5] = b_na;

        //initialize bond angles based on known constants
        b_ang[0] = ang_nac;
        b_ang[1] = ang_acn;
        b_ang[2] = ang_cna;
        b_ang[3] = ang_nac;
        b_ang[4] = ang_acn;
        b_ang[5] = ang_cna;
        b_ang[6] = ang_nac;

        //initialize peptide torsion angles based on known constants
        t_ang[0] = Math.PI;
        t_ang[1] = Math.PI;

        loopClosure.initializeLoopClosure(b_len, b_ang, t_ang);

        calc_rmsd = true;
        //write_out_pdb = true;
        //write_out_pdb = 1;
        //strcpy(in_pdb, "data/135l.pdb")

        //n1 and n2 determine the starting and ending residues of the loop in
        //      the C version of loopClosure-will likely become a user
        //      chosen critera.
        n1 = 2;
        n2 = 2;

        for (n0 = n1; n0 <= n2; n0++) {
            //choose starting and ending residues of the loop
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
                    char chainID = molAss.getBackBoneAtoms().get(i).getChainID();

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

                        //get chain ID
                        chain_n[ir] = chainID;

                       //print statements for debugging purposes
                        //System.out.println("r_n[ir][0]: "+r_n[ir][0]);
                        //System.out.println("r_n[ir][1]: "+r_n[ir][1]);
                        //System.out.println("r_n[ir][2]: "+r_n[ir][2]+"\n\n");
                    } else if (atmname.contentEquals(ca)) {
                        //Backbone alpha carbon coordinates are stored in r_a[]
                        r_a[ir][0] = coordinateArray[0];
                        r_a[ir][1] = coordinateArray[1];
                        r_a[ir][2] = coordinateArray[2];

                        //set res_name array
                        res_name[ir] = molAss.getBackBoneAtoms().get(i).getResidueName();

                        //get chain ID
                        chain_a[ir] = chainID;

                       //print statements for debugging purposes
                        //System.out.println("r_a[ir][0]: "+r_a[ir][0]);
                        //System.out.println("r_a[ir][1]: "+r_a[ir][1]);
                        //System.out.println("r_a[ir][2]: "+r_a[ir][2]+"\n\n");
                    } else if (atmname.contentEquals(c)) {
                        //Backbone carbon coordinates are stored in r_c[]
                        r_c[ir][0] = coordinateArray[0];
                        r_c[ir][1] = coordinateArray[1];
                        r_c[ir][2] = coordinateArray[2];

                        chain_c[ir] = chainID;

                       //print statements for debugging purposes
                        //System.out.println("r_c[ir][0]: "+r_c[ir][0]);
                        //System.out.println("r_c[ir][1]: "+r_c[ir][1]);
                        //System.out.println("r_c[ir][2]: "+r_c[ir][2]+"\n\n");
                    }
                    //else is only used during debugging
                    //else
                    //{
                    //System.out.println("The atom was NOT C, CA, or N.\n");
                    //}
                    i++;
                }
            }

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    //initialize the r0 arrays
                    r0_n[i][j] = r_n[i + 1][j];
                    r0_a[i][j] = r_a[i + 1][j];
                    r0_c[i][j] = r_c[i + 1][j];

                    //print statements for debugging purposes
                    //System.out.println("The r0 array " +r0_n[i][j] +" " +r0_a[i][j]+" "+r0_c[i][j] );
                }
            }

            //Method that solves the 16th degree, 3 peptide polynomial
            loopClosure.solve3PepPoly(r_n[1], r_a[1], r_a[3], r_c[3], r_soln_n, r_soln_a, r_soln_c, n_soln);

            StringBuilder sb = new StringBuilder();
            sb.append(String.format("Starting res.:             %d\n", n0));
            sb.append(String.format("Ending res.:               %d\n", (n0 + 2)));
            sb.append(String.format("No. of solutions:          %d\n", n_soln[0]));
            logger.info(sb.toString());

            //print statements for debugging purposes
            //System.out.println("Final Coordinates:\n ");
            //System.out.println("r_n: "+Arrays.deepToString(r_n)+"\n");
            //System.out.println("r_a: "+Arrays.deepToString(r_a)+"\n");
            //System.out.println("r_c: "+Arrays.deepToString(r_c)+"\n");
            //System.out.println("Calc_rmsd: "+calc_rmsd+"\n");
            //System.out.println("n_soln: "+n_soln[0]+"\n");
            if (calc_rmsd) {
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
                            dr[j] = r_soln_n[k][i][j] - r0_n[i][j];
                        }
                        sum += VectorMath.dot(dr, dr);
                        for (j = 0; j < 3; j++) {
                            dr[j] = r_soln_a[k][i][j] - r0_a[i][j];
                        }
                        sum += VectorMath.dot(dr, dr);
                        for (j = 0; j < 3; j++) {
                            dr[j] = r_soln_c[k][i][j] - r0_c[i][j];
                        }
                        sum += VectorMath.dot(dr, dr);
                    }
                    rmsd = Math.sqrt(sum / 9.0e0);

                    StringBuilder string = new StringBuilder();
                    string.append(String.format("Rmsd for solution #" + (k + 1) + " is " + rmsd + "\n"));
                    logger.info(string.toString());

                    /* if(write_out_pdb)
                     {
                     sprintf(out_pdb, "%s_%d.pdb", out_prefix, k+1);
                     */
                    counter++;
                       //if(writeFile)
                    //{
                    File fileName = sturmMethod.writePDBBackbone(out_pdb, res_name, r_n, r_a, r_c, n0 - 1, n0 + 3, chain_n, chain_a, chain_c, molAss, counter, writeFile);
                    StringBuilder string1 = new StringBuilder();
                    string.append(String.format("Recording the solution #" + (k + 1) + " in " + fileName + ".\n"));
                       //}
                    //write_pdb_backbone(out_pdb, res_name, r_n, r_a, r_c, n0-1, n0+3);

                    /*
                     writePDBBackbone(out_pdb, res_name, r_n, r_a, r_c, n0-1, n0+3);
                     }
                     */
                    /*for(res_no=stt_res;res_no<=end_res;res_no++)
                     {
                     ir = res_no - stt_res;
                     k++;
                     System.out.println("ATOM " +k+" N "+res_no+" "+ r_n[ir][0]+" "+r_n[ir][1]+" "+r_n[ir][2]);
                     k++;
                     System.out.println("ATOM " +k+" CA "+res_no+" "+ r_a[ir][0]+" "+r_a[ir][1]+" "+r_a[ir][2]);
                     k++;
                     System.out.println("ATOM " +k+" C "+res_no+" "+ r_c[ir][0]+" "+r_c[ir][1]+" "+r_c[ir][2]);
                     }
                     System.out.println("\n\n\n");*/
                }
            }
        }
    }

    // Only for JUnit testing.
    public double[][] getr_n() {
        return r_n;
    }

    public double[][] getr_a() {
        return r_a;
    }

    public double[][] getr_c() {
        return r_c;
    }
}
