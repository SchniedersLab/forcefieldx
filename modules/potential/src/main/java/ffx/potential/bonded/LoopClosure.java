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

import ffx.numerics.VectorMath;
import static ffx.numerics.VectorMath.dot;
import java.util.Arrays;
import java.util.logging.Logger;
import org.apache.commons.math3.util.FastMath;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;
//import java.util.logging.Level; may be needed

/**
 *
 * @author mrtollefson
 */
public class LoopClosure {

    /**
     *
     * @param r_n1
     * @param r_a1
     * @param r_a3
     * @param r_c3
     * @param r_soln_n
     * @param r_soln_a
     * @param r_soln_c
     * @param n_soln
     */

    private static final Logger logger = Logger.getLogger(LoopClosure.class.getName());
    private static final SturmMethod sturmMethod = new SturmMethod();

    //The initial C files from the Dill group utilize global variables. 
    //These variables should eventually be moved. 
    private int max_soln = 16;
    private int[] deg_pol = new int[1];
    private int print_level = 1;
    private double[] len0 = new double[6];
    private double[] b_ang0 = new double[7];
    private double[] t_ang0 = new double[2];
    private double aa13_min_sqr;
    private double aa13_max_sqr;
    private double[][] delta = new double[4][1];
    private double[][] xi = new double[3][1];
    private double[][] eta = new double[3][1];
    private double[] alpha = new double[3];
    private double[] theta = new double[3];
    private double[] cos_alpha = new double[3];
    private double[] sin_alpha = new double[3];
    private double[] cos_theta = new double[3];
    private double[] sin_theta = new double[3];
    private double[] cos_delta = new double[4];
    private double[] sin_delta = new double[4];
    private double[] cos_xi = new double[3];
    private double[] cos_eta = new double[3];
    private double[] sin_xi = new double[3];
    private double[] sin_eta = new double[3];
    private double[] r_a1a3 = new double[3];
    private double[] r_a1n1 = new double[3];
    private double[] r_a3c3 = new double[3];
    private double[] b_a1a3 = new double[3];
    private double[] b_a1n1 = new double[3];
    private double[] b_a3c3 = new double[3];
    private double[] len_na = new double[3];
    private double[] len_ac = new double[3];
    private double[] len_aa = new double[3];
    private double[][] C0 = new double[3][3];
    private double[][] C1 = new double[3][3];
    private double[][] C2 = new double[3][3];
    private double[][] Q = new double[5][17];
    private double[][] R = new double[3][17];
    boolean bool = true;

    public LoopClosure() 
    {
        deg_pol[0] = 16;
    }

    public void solve_3pep_poly(double[] r_n1, double[] r_a1, double[] r_a3, double[] r_c3, double[][][] r_soln_n,
            double[][][] r_soln_a, double[][][] r_soln_c, int[] n_soln) 
    {
        double[] poly_coeff = new double[deg_pol[0] + 1];
        double[] roots = new double[max_soln];
        int i;

        get_input_angles(n_soln, r_n1, r_a1, r_a3, r_c3);

        //System.out.println("n_soln: "+n_soln[0]);
        //System.out.println("r_n1: "+r_n1[0]+ " " + r_n1[1]+" "+r_n1[2]+" "+r_n1[3]+" "+r_n1[4]);
        //System.out.println("r_a1: "+r_a1[0]+ " " + r_a1[1]+" "+r_a1[2]+" "+r_a1[3]+" "+r_a1[4]);
        //System.out.println("r_a3: "+r_a3[0]+ " " + r_a3[1]+" "+r_a3[2]+" "+r_a3[3]+" "+r_a3[4]);
        //System.out.println("r_c3: "+r_c3[0]+ " " + r_c3[1]+" "+r_c3[2]+" "+r_c3[3]+" "+r_c3[4]);
        if (n_soln[0] == 0) {
            return;
        }
        
        get_poly_coeff(poly_coeff);
        //SEARCH2!!!!
        
        //System.out.println("Roots: "+Arrays.toString(roots));
        //System.out.println("\nN_soln: "+Arrays.toString(n_soln));
        //System.out.println("\ndeg_pol: "+Arrays.toString(deg_pol)+"\n"); 
        
        sturmMethod.solve_sturm(deg_pol, n_soln, poly_coeff, roots);
        
        
        if (n_soln[0] == 0) {
            System.out.println("Loop failed");
            //logger.severe("Loop building failed.");
        }

        coord_from_poly_roots(n_soln, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c);

    }

    //Returns the sign (positive or negative) of a variable.
    public double sign(double a, double b) {
        if (b >= 0.0) {
            return FastMath.abs(a);
        } else {
            return -FastMath.abs(a);
        }
    }

    public void initialize_loop_closure(double[] b_len, double[] b_ang, double[] t_ang) {
        double len1, len2, a_min, a_max;
        double[] axis1;
        double[] axis = new double[3];
        double[] rr_a1 = new double[3];
        double[] rr_c1 = new double[3];
        double[] rr_n2 = new double[3];
        double[] rr_a2 = new double[3];
        double[] rr_n2a2_ref = new double[3];
        double[] rr_c1a1 = new double[3];
        double[] rr_a1a2 = new double[3];
        double[] dr = new double[3];
        double[] bb_c1a1 = new double[3];
        double[] bb_a1a2 = new double[3];
        double[] bb_a2n2 = new double[3];
        double[] p = new double[4];
        double[][] Us = new double[3][3];
        double[] mulpro = new double[3];
        double[] tmp_val = new double[3];
        double[][] tol_secant = new double[1][1];
        int[][] max_iter_sturm = new int[1][1];
        int[][] max_iter_secant = new int[1][1];
        int i, j;

        //known constants
        tol_secant[0][0] = 1.0e-15;
        max_iter_sturm[0][0] = 100;
        max_iter_secant[0][0] = 20;

        //sets termination criteria for polynomial solver
        sturmMethod.initialize_sturm(tol_secant, max_iter_sturm, max_iter_secant);

        //initializing initial length, angle, and torsion arrays
        for (i = 0; i < 6; i++) {
            len0[i] = b_len[i];
        }

        for (i = 0; i < 7; i++) {
            b_ang0[i] = b_ang[i];
        }

        for (i = 0; i < 2; i++) {
            t_ang0[i] = t_ang[i];
        }

        for (i = 0; i < 3; i++) {
            rr_c1[i] = 0.0;
        }

        //initializing axis
        axis[0] = 1.0;
        axis[1] = 0.0;
        axis[2] = 0.0;

        for (i = 0; i < 2; i++) 
        {
            //iniitalizing rr_a1 array
            rr_a1[0] = cos(b_ang0[3 * i + 1]) * len0[3 * i];
            rr_a1[1] = sin(b_ang0[3 * i + 1]) * len0[3 * i];
            rr_a1[2] = 0.0e0;

            //initializing rr_n2 array
            rr_n2[0] = len0[3 * i + 1];
            rr_n2[1] = 0.0e0;
            rr_n2[2] = 0.0e0;

            //initializing rr_c1a1 array
            for (j = 0; j < 3; j++) 
            {
                rr_c1a1[j] = rr_a1[j] - rr_c1[j];
            }

            //initializing rr_n2a2_ref array
            rr_n2a2_ref[0] = -cos(b_ang0[3 * i + 2]) * len0[3 * i + 2];
            rr_n2a2_ref[1] = sin(b_ang0[3 * i + 2]) * len0[3 * i + 2];
            rr_n2a2_ref[2] = 0.0e0;

            //print statements for debugging purposes
            //System.out.println("rr_a1: "+rr_a1[0]+" "+rr_a1[1]+ " "+rr_a1[2]+ "\n");
            //System.out.println("rr_n2: "+rr_n2[0]+ " "+rr_n2[1]+" "+rr_n2[2]+ "\n");
            //System.out.println("rr_n2a2_ref: "+rr_n2a2_ref[0]+ " "+rr_n2a2_ref[1]+" "+rr_n2a2_ref[2]+ "\n");
            
            //quaternion is the quotient of two vectors in 3D space
            quaternion(axis, t_ang0[i] * 0.25e0, p);
            
            //means of representing a rotation of an axis about an origin
            //      with no angular specification
            rotation_matrix(p, Us);
            
            //basic matrix multiplication
            matmul(Us, rr_n2a2_ref, mulpro);

            //print statements for debugging purposes
            //System.out.println("Us: " +Us[0][0] +" "+Us[1][0] +" "+Us[2][0] +"\n");
            //System.out.println("Us: " +Us[0][1] +" "+Us[1][1] +" "+Us[2][1] +"\n");
            //System.out.println("Us: " +Us[0][2] +" "+Us[1][2] +" "+Us[2][2] +"\n");
            //System.out.println("rr_n2a2_ref: "  +rr_n2a2_ref[0] +" "+rr_n2a2_ref[1] + " " + rr_n2a2_ref[2] + "\n");
            //System.out.println("mulpro: "  +mulpro[0] +" "+mulpro[1] + " " + mulpro[2] + "\n");
            
            for (j = 0; j < 3; j++) 
            {
                rr_a2[j] = mulpro[j] + rr_n2[j];
                rr_a1a2[j] = rr_a2[j] - rr_a1[j];
                dr[j] = rr_a1a2[j];

                //print statements for debugging purposes
                //System.out.println("rr_a2: "+rr_a2[j]+"\n");
                //System.out.println("rr_a1a2: "+rr_a1a2[j]+"\n");
                //System.out.println("dr: "+dr[j]+"\n");
            }

            len2 = dot(dr, dr);
            len1 = sqrt(len2);
            len_aa[i + 1] = len1;
            
            //print statements for debugging purposes
            //System.out.println("len 2: "+len2+"\n");
            //System.out.println("len 1: "+len1+"\n");
            //System.out.println("len_aa: "+len_aa[i+1]+"\n");

            for (j = 0; j < 3; j++) {
                bb_c1a1[j] = rr_c1a1[j] / len0[3 * i];
                bb_a1a2[j] = rr_a1a2[j] / len1;
                bb_a2n2[j] = (rr_n2[j] - rr_a2[j]) / len0[3 * i + 2];

                //print statements for debugging purposes
                //System.out.println("bb_c1a1: "+bb_c1a1[j]+"\n");
                //System.out.println("bb_a1a2: "+bb_a1a2[j]+"\n");
                //System.out.println("bb_a2n2: "+bb_a2n2[j]+"\n");
            }

            for (j = 0; j < 3; j++) 
            {
                tmp_val[j] = -bb_a1a2[j];
                //System.out.println("tmp_val: "+tmp_val[j]+"\n");
            }

            //System.out.println("XI: "+xi[i+1]+"\n");
            calc_bnd_ang(tmp_val, bb_a2n2, xi[i + 1]);

            for (j = 0; j < 3; j++) 
            {
                tmp_val[j] = -bb_c1a1[j];
                //System.out.println("tmp_val: "+tmp_val[j]+"\n");
            }

            calc_bnd_ang(bb_a1a2, tmp_val, eta[i]);

            calc_dih_ang(bb_c1a1, bb_a1a2, bb_a2n2, delta[i + 1]);

            delta[i + 1][0] = FastMath.PI - delta[i + 1][0];
        }

        a_min = b_ang[3] - (xi[1][0] + eta[1][0]);
        a_max = FastMath.min((b_ang[3] + (xi[1][0] + eta[1][0])), FastMath.PI);
        aa13_min_sqr = FastMath.pow(len_aa[1], 2) + FastMath.pow(len_aa[2], 2) - 2.0e0 * len_aa[1] * len_aa[2] * cos(a_min);
        aa13_max_sqr = FastMath.pow(len_aa[1], 2) + FastMath.pow(len_aa[2], 2) - 2.0e0 * len_aa[1] * len_aa[2] * cos(a_max);

        //System.out.println("delta: "+delta[i+1][0]+"\n");
        //System.out.println("a_min: "+a_min+"\n");
        //System.out.println("a_max: "+a_max+"\n");
        //System.out.println("aa13_min_sqr: "+aa13_min_sqr+"\n");
        //System.out.println("aa13_max_sqr: "+aa13_max_sqr+"\n");
    }

    //basic 2D matrix multiplication
    //this function multiplies matrices ma and mb, storing the answer in mc.
    public void matmul(double[][] ma, double[] mb, double[] mc) 
    {
        int i, j;

        for (i = 0; i < 3; i++) 
        {
            mc[i] = 0.0;
            for (j = 0; j < 3; j++) 
            {
                mc[i] += ma[i][j] * mb[j];
            }
        }
    }

    public void get_input_angles(int[] n_soln, double[] r_n1, double[] r_a1, double[] r_a3, double[] r_c3) {
        double dr_sqr;
        double[] tmp_val = new double[3];
        int i;
        char[] cone_type = new char[2];

        n_soln[0] = max_soln;

        for (i = 0; i < 3; i++) {
            r_a1a3[i] = r_a3[i] - r_a1[i];
        }

        dr_sqr = dot(r_a1a3, r_a1a3);
        len_aa[0] = sqrt(dr_sqr);

        if ((dr_sqr < aa13_min_sqr) || (dr_sqr > aa13_max_sqr)) {
            n_soln[0] = 0;
            return;
        }

        for (i = 0; i < 3; i++) {
            r_a1n1[i] = r_n1[i] - r_a1[i];
        }

        len_na[0] = sqrt(dot(r_a1n1, r_a1n1));
        len_na[1] = len0[2];
        len_na[2] = len0[5];

        for (i = 0; i < 3; i++) {
            r_a3c3[i] = r_c3[i] - r_a3[i];
        }

        len_ac[0] = len0[0];
        len_ac[1] = len0[3];
        len_ac[2] = sqrt(dot(r_a3c3, r_a3c3));

        for (i = 0; i < 3; i++) {
            b_a1n1[i] = r_a1n1[i] / len_na[0];
            b_a3c3[i] = r_a3c3[i] / len_ac[2];
            b_a1a3[i] = r_a1a3[i] / len_aa[0];
        }

        for (i = 0; i < 3; i++) {
            tmp_val[i] = -b_a1n1[i];
        }

        calc_dih_ang(tmp_val, b_a1a3, b_a3c3, delta[3]);

        delta[0][0] = delta[3][0];

        for (i = 0; i < 3; i++) {
            tmp_val[i] = -b_a1a3[i];
        }

        calc_bnd_ang(tmp_val, b_a1n1, xi[0]);

        calc_bnd_ang(b_a1a3, b_a3c3, eta[2]);

        for (i = 0; i < 3; i++) {
            cos_delta[i + 1] = cos(delta[i + 1][0]);
            sin_delta[i + 1] = sin(delta[i + 1][0]);
            cos_xi[i] = cos(xi[i][0]);
            sin_xi[i] = sin(xi[i][0]);
            sin_xi[i] = sin(xi[i][0]);
            cos_eta[i] = cos(eta[i][0]);
            cos_eta[i] = cos(eta[i][0]);
            sin_eta[i] = sin(eta[i][0]);
            sin_eta[i] = sin(eta[i][0]);
        }

        cos_delta[0] = cos_delta[3];
        sin_delta[0] = sin_delta[3];
        theta[0] = b_ang0[0];
        theta[1] = b_ang0[3];
        theta[2] = b_ang0[6];

        for (i = 0; i < 3; i++) {
            cos_theta[i] = cos(theta[i]);
        }

        cos_alpha[0] = -((len_aa[0] * len_aa[0]) + (len_aa[1] * len_aa[1]) - (len_aa[2] * len_aa[2])) / (2.0e0 * len_aa[0] * len_aa[1]);
        alpha[0] = acos(cos_alpha[0]);
        sin_alpha[0] = sin(alpha[0]);
        cos_alpha[1] = ((len_aa[1] * len_aa[1]) + (len_aa[2] * len_aa[2]) - (len_aa[0] * len_aa[0])) / (2.0e0 * len_aa[1] * len_aa[2]);
        alpha[1] = acos(cos_alpha[1]);
        sin_alpha[1] = sin(alpha[1]);
        alpha[2] = FastMath.PI - alpha[0] + alpha[1];
        cos_alpha[2] = cos(alpha[2]);
        sin_alpha[2] = sin(alpha[2]);

        //if (print_level > 0)
        //{
        logger.info(String.format("xi = %9.4f%9.4f%9.4f\\n", toDegrees(xi[0][0]), toDegrees(xi[1][0]), toDegrees(xi[2][0])));
        logger.info(String.format("eta = %9.4f%9.4f%9.4f\n", toDegrees(eta[0][0]), toDegrees(eta[1][0]), toDegrees(eta[2][0])));
        logger.info(String.format("delta = %9.4f%9.4f%9.4f\n", toDegrees(delta[1][0]), toDegrees(delta[2][0]), toDegrees(delta[3][0])));
        logger.info(String.format("theta = %9.4f%9.4f%9.4f\n", toDegrees(theta[0]), toDegrees(theta[1]), toDegrees(theta[2])));
        logger.info(String.format("alpha = %9.4f%9.4f%9.4f\n", toDegrees(alpha[0]), toDegrees(alpha[1]), toDegrees(alpha[2])));
       // }

        for (i = 0; i < 3; i++) {
            test_two_cone_existence_soln(theta[i], xi[i][0], eta[i][0], alpha[i], n_soln, cone_type);
            if (n_soln[0] == 0) {
                return;
            }
        }

    }

    public void test_two_cone_existence_soln(double tt, double kx, double et, double ap, int[] n_soln, char[] cone_type) {
        double at;
        double ex;
        double abs_at;
        double ap1;
        double kx1;
        double et1;
        double cos_tx1;
        double cos_tx2;
        double cos_te1;
        double cos_te2;
        double cos_ea1;
        double cos_ea2;
        double cos_xa1;
        double cos_xa2;
        int s1;
        int s2;
        int t1;
        int t2;
        boolean complicated = false;

        n_soln[0] = max_soln;
        ap1 = ap;
        kx1 = kx;
        et1 = et;

        at = ap1 - tt;
        ex = kx1 + et1;
        abs_at = FastMath.abs(at);

        if (abs_at > ex) {
            n_soln[0] = 0;
            return;
        }

        if (complicated) {
            cos_tx1 = cos(tt + kx1);
            cos_tx2 = cos(tt - kx1);
            cos_te1 = cos(tt + et1);
            cos_te2 = cos(tt - et1);
            cos_ea1 = cos(et1 + ap1);
            cos_ea2 = cos(et1 - ap1);
            cos_xa1 = cos(kx1 + ap1);
            cos_xa2 = cos(kx1 - ap1);

            s1 = 0;
            s2 = 0;
            t1 = 0;
            t2 = 0;

            if ((cos_te1 - cos_xa2) * (cos_te1 - cos_xa1) <= 0.0e0) {
                s1 = 0;
            }

            if ((cos_te2 - cos_xa2) * (cos_te2 - cos_xa1) <= 0.0e0) {
                s2 = 0;
            }

            if ((cos_tx1 - cos_ea2) * (cos_tx1 - cos_ea1) <= 0.0e0) {
                t1 = 0;
            }

            if ((cos_tx2 - cos_ea2) * (cos_tx2 - cos_ea1) <= 0.0e0) {
                t2 = 0;
            }
        }
    }

    public void get_poly_coeff(double[] poly_coeff) //poly_coeff[deg_pol+1]
    {
        int i;
        int j;
        double A0, A1, A2, A3, A4, A21, A22, A31, A32, A41, A42;
        double[] B0 = new double[3];
        double[] B1 = new double[3];
        double[] B2 = new double[3];
        double[] B3 = new double[3];
        double[] B4 = new double[3];
        double[] B5 = new double[3];
        double[] B6 = new double[3];
        double[] B7 = new double[3];
        double[] B8 = new double[3];
        double[][] u11 = new double[5][5];
        double[][] u12 = new double[5][5];
        double[][] u13 = new double[5][5];
        double[][] u31 = new double[5][5];
        double[][] u32 = new double[5][5];
        double[][] u33 = new double[5][5];
        double[][] um1 = new double[5][5];
        double[][] um2 = new double[5][5];
        double[][] um3 = new double[5][5];
        double[][] um4 = new double[5][5];
        double[][] um5 = new double[5][5];
        double[][] um6 = new double[5][5];
        double[][] q_tmp = new double[5][5];
        int[] p1 = new int[2];
        int[] p3 = new int[2];
        int[] p_um1 = new int[2];
        int[] p_um2 = new int[2];
        int[] p_um3 = new int[2];
        int[] p_um4 = new int[2];
        int[] p_um5 = new int[2];
        int[] p_um6 = new int[2];
        int[] p_Q = new int[2];
        int p2, p4;
        int[] p_f1 = new int[1];
        int[] p_f2 = new int[1];
        int[] p_f3 = new int[1];
        int[] p_f4 = new int[1];
        int[] p_f5 = new int[1];
        int[] p_f6 = new int[1];
        int[] p_f7 = new int[1];
        int[] p_f8 = new int[1];
        int[] p_f9 = new int[1];
        int[] p_f10 = new int[1];
        int[] p_f11 = new int[1];
        int[] p_f12 = new int[1];
        int[] p_f13 = new int[1];
        int[] p_f14 = new int[1];
        int[] p_f15 = new int[1];
        int[] p_f16 = new int[1];
        int[] p_f17 = new int[1];
        int[] p_f18 = new int[1];
        int[] p_f19 = new int[1];
        int[] p_f20 = new int[1];
        int[] p_f21 = new int[1];
        int[] p_f22 = new int[1];
        int[] p_f23 = new int[1];
        int[] p_f24 = new int[1];
        int[] p_f25 = new int[1];
        int[] p_f26 = new int[1];
        int[] p_final = new int[1];
        double[] f1 = new double[17];
        double[] f2 = new double[17];
        double[] f3 = new double[17];
        double[] f4 = new double[17];
        double[] f5 = new double[17];
        double[] f6 = new double[17];
        double[] f7 = new double[17];
        double[] f8 = new double[17];
        double[] f9 = new double[17];
        double[] f10 = new double[17];
        double[] f11 = new double[17];
        double[] f12 = new double[17];
        double[] f13 = new double[17];
        double[] f14 = new double[17];
        double[] f15 = new double[17];
        double[] f16 = new double[17];
        double[] f17 = new double[17];
        double[] f18 = new double[17];
        double[] f19 = new double[17];
        double[] f20 = new double[17];
        double[] f21 = new double[17];
        double[] f22 = new double[17];
        double[] f23 = new double[17];
        double[] f24 = new double[17];
        double[] f25 = new double[17];
        double[] f26 = new double[17];

        for (i = 0; i < 3; i++) {
            A0 = cos_alpha[i] * cos_xi[i] * cos_eta[i] - cos_theta[i];
            A1 = -sin_alpha[i] * cos_xi[i] * sin_eta[i];
            A2 = sin_alpha[i] * sin_xi[i] * cos_eta[i];
            A3 = sin_xi[i] * sin_eta[i];
            A4 = A3 * cos_alpha[i];
            j = i;
            A21 = A2 * cos_delta[j];
            A22 = A2 * sin_delta[j];
            A31 = A3 * cos_delta[j];
            A32 = A3 * sin_delta[j];
            A41 = A4 * cos_delta[j];
            A42 = A4 * sin_delta[j];
            B0[i] = A0 + A22 + A31;
            B1[i] = 2.0 * (A1 + A42);
            B2[i] = 2.0 * (A32 - A21);
            B3[i] = -4.0 * A41;
            B4[i] = A0 + A22 - A31;
            B5[i] = A0 - A22 - A31;
            B6[i] = -2.0 * (A21 + A32);
            B7[i] = 2.0 * (A1 - A42);
            B8[i] = A0 - A22 + A31;

            //print statements for debugging purposes
            //System.out.print("A0: "+A0+"\n");
            //System.out.println("A1: "+A1+"\n");
            //System.out.println("A2: "+A2+"\n");
            //System.out.println("A3: "+A3+"\n");
            //System.out.println("A4: "+A4+"\n");
            //System.out.println("A21: "+A21+"\n");
            //System.out.println("A22: "+A22+"\n");
            //System.out.println("A31: "+A31+"\n");
            //System.out.println("A32: "+A32+"\n");
            //System.out.println("A41: "+A41+"\n");
            //System.out.println("A42: "+A42+"\n");    
        }
           //print statements for debugging purposes
           //System.out.println("B0: "+Arrays.toString(B0)+"\n\n");
           //System.out.println("B1: "+Arrays.toString(B1)+"\n\n");
           //System.out.println("B2: "+Arrays.toString(B2)+"\n\n");
           //System.out.println("B3: "+Arrays.toString(B3)+"\n\n");
           //System.out.println("B4: "+Arrays.toString(B4)+"\n\n");
           //System.out.println("B5: "+Arrays.toString(B5)+"\n\n");
           //System.out.println("B6: "+Arrays.toString(B6)+"\n\n");
           //System.out.println("B7: "+Arrays.toString(B7)+"\n\n");
           //System.out.println("B8: "+Arrays.toString(B8)+"\n\n");

        i = 0;

        C0[i][0] = B0[i];
        C0[i][1] = B2[i];
        C0[i][2] = B5[i];

        C1[i][0] = B1[i];
        C1[i][1] = B3[i];
        C1[i][2] = B7[i];

        C2[i][0] = B4[i];
        C2[i][1] = B6[i];
        C2[i][2] = B8[i];

        for (i = 1; i < 3; i++) {
            C0[i][0] = B0[i];
            C0[i][1] = B1[i];
            C0[i][2] = B4[i];

            C1[i][0] = B2[i];
            C1[i][1] = B3[i];
            C1[i][2] = B6[i];

            C2[i][0] = B5[i];
            C2[i][1] = B7[i];
            C2[i][2] = B8[i];
        }
        
        //print statements for debugging purposes
        //System.out.println("C0: "+Arrays.deepToString(C0));
        //System.out.println("C1: "+Arrays.deepToString(C1));
        //System.out.println("C2: "+Arrays.deepToString(C2));

        for (i = 0; i < 3; i++) {
            u11[0][i] = C0[0][i];
            u12[0][i] = C1[0][i];
            u13[0][i] = C2[0][i];
            u31[i][0] = C0[1][i];
            u32[i][0] = C1[1][i];
            u33[i][0] = C2[1][i];
            //System.out.println("U31: "+u31[i][0]+"\n");
        }

        p1[0] = 2;
        p1[1] = 0;
        p3[0] = 0;
        p3[1] = 2;
        
        //print statements for debugging purposes
        //System.out.println("U31: "+Arrays.deepToString(u31)+"\n\n\n");
        //System.out.println("p3: "+Arrays.toString(p3)+"\n\n\n");
        
        poly_mul_sub2(u32, u32, u31, u33, p3, p3, p3, p3, um1, p_um1);
        
        //print statements for debugging purposes
        //System.out.println("p3: "+Arrays.toString(p3)+"\n\n\n");
        //System.out.println("U31: "+Arrays.deepToString(u31)+"\n\n\n");
        //System.out.println("U33: "+Arrays.deepToString(u33)+"\n\n\n");
        
        poly_mul_sub2(u12, u32, u11, u33, p1, p3, p1, p3, um2, p_um2);
        poly_mul_sub2(u12, u33, u13, u32, p1, p3, p1, p3, um3, p_um3);
        poly_mul_sub2(u11, u33, u31, u13, p1, p3, p3, p1, um4, p_um4);
        poly_mul_sub2(u13, um1, u33, um2, p1, p_um1, p3, p_um2, um5, p_um5);
        poly_mul_sub2(u13, um4, u12, um3, p1, p_um4, p1, p_um3, um6, p_um6);
        poly_mul_sub2(u11, um5, u31, um6, p1, p_um5, p3, p_um6, q_tmp, p_Q);

        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Q[i][j] = q_tmp[i][j];
            }
        }

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 17; j++) {
                R[i][j] = 0.0;
            }
        }

        for (i = 0; i < 3; i++) {
            R[0][i] = C0[2][i];
            R[1][i] = C1[2][i];
            R[2][i] = C2[2][i];
            
            //print statements for debugging purposes
            //System.out.println("R[0][i]: "+R[0][i]);
            //System.out.println("R[1][i]: "+R[1][i]);
            //System.out.println("R[2][i]: "+R[2][i]);
        }

        p2 = 2;
        p4 = 4;

        poly_mul_sub1(R[1], R[1], R[0], R[2], p2, p2, p2, p2, f1, p_f1);
        poly_mul1(R[1], R[2], p2, p2, f2, p_f2);
        poly_mul_sub1(R[1], f1, R[0], f2, p2, p_f1[0], p2, p_f2[0], f3, p_f3);
        poly_mul1(R[2], f1, p2, p_f1[0], f4, p_f4);
        poly_mul_sub1(R[1], f3, R[0], f4, p2, p_f3[0], p2, p_f4[0], f5, p_f5);
        poly_mul_sub1(Q[1], R[1], Q[0], R[2], p4, p2, p4, p2, f6, p_f6);
        poly_mul_sub1(Q[2], f1, R[2], f6, p4, p_f1[0], p2, p_f6[0], f7, p_f7);
        poly_mul_sub1(Q[3], f3, R[2], f7, p4, p_f3[0], p2, p_f7[0], f8, p_f8);
        poly_mul_sub1(Q[4], f5, R[2], f8, p4, p_f5[0], p2, p_f8[0], f9, p_f9);
        poly_mul_sub1(Q[3], R[1], Q[4], R[0], p4, p2, p4, p2, f10, p_f10);
        poly_mul_sub1(Q[2], f1, R[0], f10, p4, p_f1[0], p2, p_f10[0], f11, p_f11);
        poly_mul_sub1(Q[1], f3, R[0], f11, p4, p_f3[0], p2, p_f11[0], f12, p_f12);
        poly_mul_sub1(Q[2], R[1], Q[1], R[2], p4, p2, p4, p2, f13, p_f13);
        poly_mul_sub1(Q[3], f1, R[2], f13, p4, p_f1[0], p2, p_f13[0], f14, p_f14);
        poly_mul_sub1(Q[3], R[1], Q[2], R[2], p4, p2, p4, p2, f15, p_f15);
        poly_mul_sub1(Q[4], f1, R[2], f15, p4, p_f1[0], p2, p_f15[0], f16, p_f16);
        poly_mul_sub1(Q[1], f14, Q[0], f16, p4, p_f14[0], p4, p_f16[0], f17, p_f17);
        poly_mul_sub1(Q[2], R[2], Q[3], R[1], p4, p2, p4, p2, f18, p_f18);
        poly_mul_sub1(Q[1], R[2], Q[3], R[0], p4, p2, p4, p2, f19, p_f19);
        poly_mul_sub1(Q[3], f19, Q[2], f18, p4, p_f19[0], p4, p_f18[0], f20, p_f20);
        poly_mul_sub1(Q[1], R[1], Q[2], R[0], p4, p2, p4, p2, f21, p_f21);
        poly_mul1(Q[4], f21, p4, p_f21[0], f22, p_f22);
        poly_sub1(f20, f22, p_f20, p_f22, f23, p_f23);
        poly_mul1(R[0], f23, p2, p_f23[0], f24, p_f24);
        poly_sub1(f17, f24, p_f17, p_f24, f25, p_f25);
        poly_mul_sub1(Q[4], f12, R[2], f25, p4, p_f12[0], p2, p_f25[0], f26, p_f26);
        poly_mul_sub1(Q[0], f9, R[0], f26, p4, p_f9[0], p2, p_f26[0], poly_coeff, p_final); 
        
        //print statements for debugging purposes
        //System.out.println("R[1]: "+Arrays.toString(R[1]));
        //System.out.println("R[2]: "+Arrays.toString(R[2]));
        //System.out.println("p2: "+p2);
        //System.out.println("f2: "+Arrays.toString(f2));
        //System.out.println("p_f2: "+Arrays.toString(p_f2));

        if (p_final[0] != 16) {
            System.out.println("Error. Degree of polynomial is not 16!\n");
            return;
        }

        if (poly_coeff[16] < 0.0e0) {
            for (i = 0; i < 17; i++) {
                poly_coeff[i] *= -1.0;
            }
        }

        if (print_level > 0) {
            System.out.println("poly_coeff\n");
            for (i = 0; i < 17; i++) {
                logger.info(String.format("%5d%15.6f\n", i, poly_coeff[i]));
            }
        }

        /* System.out.println("R: "+Arrays.deepToString(R)+"\n");
         System.out.println("Q: "+Arrays.deepToString(Q)+"\n");
         System.out.println("B0: "+Arrays.toString(B0)+"\n");
         System.out.println("B1: "+Arrays.toString(B1)+"\n");
         System.out.println("B2: "+Arrays.toString(B2)+"\n");
         System.out.println("B3: "+Arrays.toString(B3)+"\n");
         System.out.println("B4: "+Arrays.toString(B4)+"\n");
         System.out.println("B5: "+Arrays.toString(B5)+"\n");
         System.out.println("B6: "+Arrays.toString(B6)+"\n");
         System.out.println("B7: "+Arrays.toString(B7)+"\n");
         System.out.println("B8: "+Arrays.toString(B8)+"\n");
         System.out.println("um1: " + Arrays.deepToString(um1) + "\n\n\n");
         System.out.println("um2: " + Arrays.deepToString(um2) + "\n\n\n");
         System.out.println("um3: " + Arrays.deepToString(um3) + "\n\n\n");
         System.out.println("um4: " + Arrays.deepToString(um4) + "\n\n\n");
         System.out.println("um5: " + Arrays.deepToString(um5) + "\n\n\n");
         System.out.println("um6: "+Arrays.deepToString(um6)+"\n");
         System.out.println("q_tmp: "+Arrays.deepToString(q_tmp)+"\n");
         System.out.println("pQ: "+Arrays.toString(p_Q)+"\n");
         System.out.println("f1: "+Arrays.toString(f1)+"\n");
         System.out.println("f2: "+Arrays.toString(f2)+"\n");
         System.out.println("f3: "+Arrays.toString(f3)+"\n");
         System.out.println("f4: "+Arrays.toString(f4)+"\n");
         System.out.println("f5: "+Arrays.toString(f5)+"\n");
         System.out.println("f6: "+Arrays.toString(f6)+"\n");
         System.out.println("f7: "+Arrays.toString(f7)+"\n");
         System.out.println("f8: "+Arrays.toString(f8)+"\n");
         System.out.println("f9: "+Arrays.toString(f9)+"\n");
         System.out.println("f10: "+Arrays.toString(f10)+"\n");
         System.out.println("f11: "+Arrays.toString(f11)+"\n");
         System.out.println("f12: "+Arrays.toString(f12)+"\n");
         System.out.println("f13: "+Arrays.toString(f13)+"\n");
         System.out.println("f14: "+Arrays.toString(f14)+"\n");
         System.out.println("f15: "+Arrays.toString(f15)+"\n");
         System.out.println("f16: "+Arrays.toString(f16)+"\n");
         System.out.println("f17: "+Arrays.toString(f17)+"\n");
         System.out.println("f18: "+Arrays.toString(f18)+"\n");
         System.out.println("f19: "+Arrays.toString(f19)+"\n");
         System.out.println("f20: "+Arrays.toString(f20)+"\n");
         System.out.println("f21: "+Arrays.toString(f21)+"\n");
         System.out.println("f22: "+Arrays.toString(f22)+"\n");
         System.out.println("f23: "+Arrays.toString(f23)+"\n");
         System.out.println("f24: "+Arrays.toString(f24)+"\n");
         System.out.println("f25: "+Arrays.toString(f25)+"\n");
         System.out.println("f26: "+Arrays.toString(f26)+"\n"); */

    }

    public void poly_mul_sub2(double[][] u1, double[][] u2, double[][] u3, double[][] u4, int[] p1, int[] p2, int[] p3, int[] p4, double[][] u5, int[] p5) {
        double[][] d1 = new double[5][5];
        double[][] d2 = new double[5][5];
        int[] pd1 = new int[2];
        int[] pd2 = new int[2];

        poly_mul2(u1, u2, p1, p2, d1, pd1);
        poly_mul2(u3, u4, p3, p4, d2, pd2);
        poly_sub2(d1, d2, pd1, pd2, u5, p5);

        /*int i;
        int j;
        for (i=0; i<5;i++)
         {
             for (j=0;j<5;j++)
             {
                 System.out.println(u5[i][j]+"\n");
             }
             
         } */
    }

    public void poly_mul2(double[][] u1, double[][] u2, int[] p1, int[] p2, double[][] u3, int[] p3) {
        int i1, j1, i2, j2, i3, j3, p11, p12, p21, p22;
        int i, j;
        double u1ij;

        for (i = 0; i < 2; i++) {
            p3[i] = p1[i] + p2[i];
        }

        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                u3[i][j] = 0.0;
            }
        }

        p11 = p1[0];
        p12 = p1[1];
        p21 = p2[0];
        p22 = p2[1];

        for (i1 = 0; i1 <= p12; i1++) {
            for (j1 = 0; j1 <= p11; j1++) {
                u1ij = u1[i1][j1];
                //System.out.println("u1ij: "+u1ij+"\n");
                for (i2 = 0; i2 <= p22; i2++) {
                    i3 = i1 + i2;
                    for (j2 = 0; j2 <= p21; j2++) {
                        j3 = j1 + j2;
                        u3[i3][j3] = u3[i3][j3] + u1ij * u2[i2][j2];
                    }
                }
            }
        }
    }

    public void poly_sub2(double[][] u1, double[][] u2, int[] p1, int[] p2, double[][] u3, int[] p3) {
        int i, j, p11, p12, p21, p22;
        boolean i1_ok, i2_ok;

        p11 = p1[0];
        p12 = p1[1];
        p21 = p2[0];
        p22 = p2[1];
        p3[0] = FastMath.max(p11, p21);
        //System.out.println("p11, p21: "+p11+" "+p21+"\n");
        //System.out.println("p3[0]: "+p3[0]+"\n");
        p3[1] = FastMath.max(p12, p22);
        //System.out.println("p12, p22: "+p12+" "+p22+"\n");
        //System.out.println("p3[1]: "+p3[1]+"\n\n\n\n");

        //System.out.println("U2: "+Arrays.deepToString(u2));
        for (i = 0; i <= p3[1]; i++) {
            i1_ok = (i > p12);
            i2_ok = (i > p22);
           //System.out.println("Value: "+i1_ok+" "+i2_ok);

            for (j = 0; j <= p3[0]; j++) {
                if (i2_ok || (j > p21)) {
                    u3[i][j] = u1[i][j];
                } else if (i1_ok || (j > p11)) {
                    u3[i][j] = -u2[i][j];
                } else {
                    u3[i][j] = u1[i][j] - u2[i][j];
                }

            }
        }

        /*for (i=0;i<5;i++)
        {
            for (j=0;j<5;j++)
            {
                System.out.println("u3: "+u3[i][j]);
            }
         } */
    }

    public void poly_mul_sub1(double[] u1, double[] u2, double[] u3, double[] u4, int p1, int p2, int p3, int p4, double[] u5, int[] p5) {
        double[] d1 = new double[17];
        double[] d2 = new double[17];
        int[] pd1 = new int[1];
        int[] pd2 = new int[1];

        poly_mul1(u1, u2, p1, p2, d1, pd1);
        poly_mul1(u3, u4, p3, p4, d2, pd2);
        poly_sub1(d1, d2, pd1, pd2, u5, p5);
    }

    public void poly_mul1(double[] u1, double[] u2, int p1, int p2, double[] u3, int[] p3) {
        int i, i1, i2, i3;
        double u1i;

        p3[0] = p1 + p2;

        for (i = 0; i < 17; i++) {
            u3[i] = 0.0;
        }

        for (i1 = 0; i1 <= p1; i1++) {
            u1i = u1[i1];
            for (i2 = 0; i2 <= p2; i2++) {
                i3 = i1 + i2;
                u3[i3] = u3[i3] + u1i * u2[i2];
            }
        }
    }

    public void poly_sub1(double[] u1, double[] u2, int[] p1, int[] p2, double[] u3, int[] p3) {
        int i;

        p3[0] = FastMath.max(p1[0], p2[0]);
        for (i = 0; i <= p3[0]; i++) {
            if (i > p2[0]) {
                u3[i] = u1[i];
            } else if (i > p1[0]) {
                u3[i] = -u2[i];
            } else {
                u3[i] = u1[i] - u2[i];
            }
        }
    }

    public void coord_from_poly_roots(int[] n_soln, double[] roots, double[] r_n1, double[] r_a1, double[] r_a3, double[] r_c3, double[][][] r_soln_n, double[][][] r_soln_a, double[][][] r_soln_c) {
        double[] ex = new double[3];
        double[] ey = new double[3];
        double[] ez = new double[3];
        double[] b_a1a2 = new double[3];
        double[] b_a3a2 = new double[3];
        double[] r_tmp = new double[3];

        double[][] p_s = new double[3][3];
        double[][] s1 = new double[3][3];
        double[][] s2 = new double[3][3];
        double[][] p_t = new double[3][3];
        double[][] t1 = new double[3][3];
        double[][] t2 = new double[3][3];
        double[][] p_s_c = new double[3][3];
        double[][] s1_s = new double[3][3];
        double[][] s2_s = new double[3][3];
        double[][] p_t_c = new double[3][3];
        double[][] t1_s = new double[3][3];
        double[][] t2_s = new double[3][3];
        double[] angle = new double[1];
        double sig1_init;
        double[] half_tan = new double[3];
        double[] cos_tau = new double[4];
        double[] sin_tau = new double[4];
        double[] cos_sig = new double[3];
        double[] sin_sig = new double[3];
        double ht, tmp, sig1;
        double[] r_s = new double[3];
        double[] r_t = new double[3];
        double[] r0 = new double[3];
        double[][] r_n = new double[3][3];
        double[][] r_a = new double[3][3];
        double[][] r_c = new double[3][3];
        double[] p = new double[4];
        double[][] Us = new double[3][3];
        int i_soln, i, j;
        double a1c1, c1n2, n2a2, a2c2, c2n3, n3a3, a1a2, a2a3;
        double[] rr_a1c1 = new double[3];
        double[] rr_c1n2 = new double[3];
        double[] rr_n2a2 = new double[3];
        double[] rr_a2c2 = new double[3];
        double[] rr_c2n3 = new double[3];
        double[] rr_n3a3 = new double[3];
        double[] rr_a1a2 = new double[3];
        double[] rr_a2a3 = new double[3];
        double[] a3a1a2 = new double[1];
        double[] a2a3a1 = new double[1];
        double[] n1a1c1 = new double[1];
        double[] n2a2c2 = new double[1];
        double[] n3a3c3 = new double[1];
        double[] a1c1n2a2 = new double[1];
        double[] a2c2n3a3 = new double[1];
        double tmp_value;
        double[] ex_tmp = new double[3];
        double[] tmp_array = new double[3];
        double[] tmp_array1 = new double[3];
        double[] tmp_array2 = new double[3];
        double[] tmp_array3 = new double[3];
        double[] mat1 = new double[3];
        double[] mat2 = new double[3];
        double[] mat3 = new double[3];
        double[] mat4 = new double[3];
        double[] mat5 = new double[3];
        double[] mat11 = new double[3];
        double[] mat22 = new double[3];
        double[] mat33 = new double[3];
        double[] mat44 = new double[3];
        double[] mat55 = new double[3];

        if (n_soln[0] == 0) {
            return;
        }

        for (i = 0; i < 3; i++) {
            ex[i] = b_a1a3[i];
        }

        VectorMath.cross(r_a1n1, ex, ez);

        tmp_value = FastMath.sqrt(VectorMath.dot(ez, ez));

        for (i = 0; i < 3; i++) {
            ez[i] = ez[i] / tmp_value;
        }

        VectorMath.cross(ez, ex, ey);

        for (i = 0; i < 3; i++) {
            b_a1a2[i] = -cos_alpha[0] * ex[i] + sin_alpha[0] * ey[i];
            b_a3a2[i] = cos_alpha[2] * ex[i] + sin_alpha[2] * ey[i];
        }

        for (i = 0; i < 3; i++) {
            p_s[0][i] = -ex[i];
            s1[0][i] = ez[i];
            s2[0][i] = ey[i];
            p_t[0][i] = b_a1a2[i];
            t1[0][i] = ez[i];
            t2[0][i] = sin_alpha[0] * ex[i] + cos_alpha[0] * ey[i];
        }

        for (i = 0; i < 3; i++) {
            p_s[1][i] = -b_a1a2[i];
            s1[1][i] = -ez[i];
            s2[1][i] = t2[0][i];
            p_t[1][i] = -b_a3a2[i];
            t1[1][i] = -ez[i];
            t2[1][i] = sin_alpha[2] * ex[i] - cos_alpha[2] * ey[i];
        }

        for (i = 0; i < 3; i++) {
            p_s[2][i] = b_a3a2[i];
            s2[2][i] = t2[1][i];
            s1[2][i] = ez[i];
            p_t[2][i] = ex[i];
            t1[2][i] = ez[i];
            t2[2][i] = -ey[i];
        }

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                p_s_c[i][j] = p_s[i][j] * cos_xi[i];
                s1_s[i][j] = s1[i][j] * sin_xi[i];
                s2_s[i][j] = s2[i][j] * sin_xi[i];
                p_t_c[i][j] = p_t[i][j] * cos_eta[i];
                t1_s[i][j] = t1[i][j] * sin_eta[i];
                t2_s[i][j] = t2[i][j] * sin_eta[i];
            }
        }

        for (i = 0; i < 3; i++) {
            r_tmp[i] = (r_a1n1[i] / len_na[0] - p_s_c[0][i]) / sin_xi[0];
        }

        calc_bnd_ang(s1[0], r_tmp, angle);

        sig1_init = sign(angle[0], dot(r_tmp, s2[0]));

        for (i = 0; i < 3; i++) {
            r_a[0][i] = r_a1[i];
            r_a[1][i] = r_a1[i] + len_aa[1] * b_a1a2[i];
            r_a[2][i] = r_a3[i];
            r0[i] = r_a1[i];
        }

        for (i_soln = 0; i_soln < n_soln[0]; i_soln++) {
            half_tan[2] = roots[i_soln];
            half_tan[1] = calc_t2(half_tan[2]);
            half_tan[0] = calc_t1(half_tan[2], half_tan[1]);

            for (i = 1; i <= 3; i++) {
                ht = half_tan[i - 1];
                tmp = 1.0e0 + ht * ht;
                cos_tau[i] = (1.0e0 - ht * ht) / tmp;
                sin_tau[i] = 2.0e0 * ht / tmp;
            }

            cos_tau[0] = cos_tau[3];
            sin_tau[0] = sin_tau[3];

            for (i = 0; i < 3; i++) {
                cos_sig[i] = cos_delta[i] * cos_tau[i] + sin_delta[i] * sin_tau[i];
                sin_sig[i] = sin_delta[i] * cos_tau[i] - cos_delta[i] * sin_tau[i];
            }

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    r_s[j] = p_s_c[i][j] + cos_sig[i] * s1_s[i][j] + sin_sig[i] * s2_s[i][j];
                    r_t[j] = p_t_c[i][j] + cos_tau[i + 1] * t1_s[i][j] + sin_tau[i + 1] * t2_s[i][j];
                    r_n[i][j] = r_s[j] * len_na[i] + r_a[i][j];
                    r_c[i][j] = r_t[j] * len_ac[i] + r_a[i][j];
                }
            }
            //sig1 = atan2(sin_sig[0], cos_sig[0]);

            sig1 = FastMath.atan2(sin_sig[0], cos_sig[0]);
            ex_tmp[0] = -ex[0];
            ex_tmp[1] = -ex[1];
            ex_tmp[2] = -ex[2];
            tmp_value = -(sig1 - sig1_init) * 0.25;
            quaternion(ex_tmp, tmp_value, p);

            rotation_matrix(p, Us);

            for (i = 0; i < 3; i++) {
                mat11[i] = r_c[0][i] - r0[i];
                mat22[i] = r_n[1][i] - r0[i];
                mat33[i] = r_a[1][i] - r0[i];
                mat44[i] = r_c[1][i] - r0[i];
                mat55[i] = r_n[2][i] - r0[i];
            }

            matmul(Us, mat11, mat1);
            matmul(Us, mat22, mat2);
            matmul(Us, mat33, mat3);
            matmul(Us, mat44, mat4);
            matmul(Us, mat55, mat5);

            for (i = 0; i < 3; i++) {
                r_soln_n[i_soln][0][i] = r_n1[i];
                r_soln_a[i_soln][0][i] = r_a1[i];
                r_soln_c[i_soln][0][i] = mat1[i] + r0[i];
                r_soln_n[i_soln][1][i] = mat2[i] + r0[i];
                r_soln_a[i_soln][1][i] = mat3[i] + r0[i];
                r_soln_c[i_soln][1][i] = mat4[i] + r0[i];
                r_soln_n[i_soln][2][i] = mat5[i] + r0[i];
                r_soln_a[i_soln][2][i] = r_a3[i];
                r_soln_c[i_soln][2][i] = r_c3[i];
            }

            if (print_level > 0) {
                logger.info(String.format("roots: t0, t2, t1 %d\n", i_soln));
                logger.info(String.format("%15.6f %15.6f %15.6f\n", half_tan[2], half_tan[1], half_tan[0]));
            }

            for (i = 0; i < 3; i++) {
                rr_a1c1[i] = r_soln_c[i_soln][0][i] - r_soln_a[i_soln][0][i];
                rr_c1n2[i] = r_soln_n[i_soln][1][i] - r_soln_c[i_soln][0][i];
                rr_n2a2[i] = r_soln_a[i_soln][1][i] - r_soln_n[i_soln][1][i];
                rr_a2c2[i] = r_soln_c[i_soln][1][i] - r_soln_a[i_soln][1][i];
                rr_c2n3[i] = r_soln_n[i_soln][2][i] - r_soln_c[i_soln][1][i];
                rr_n3a3[i] = r_soln_a[i_soln][2][i] - r_soln_n[i_soln][2][i];
                rr_a1a2[i] = r_soln_a[i_soln][1][i] - r_soln_a[i_soln][0][i];
                rr_a2a3[i] = r_soln_a[i_soln][2][i] - r_soln_a[i_soln][1][i];
            }

            a1c1 = sqrt(dot(rr_a1c1, rr_a1c1));
            c1n2 = sqrt(dot(rr_c1n2, rr_c1n2));
            n2a2 = sqrt(dot(rr_n2a2, rr_n2a2));
            a2c2 = sqrt(dot(rr_a2c2, rr_a2c2));
            c2n3 = sqrt(dot(rr_c2n3, rr_c2n3));
            n3a3 = sqrt(dot(rr_n3a3, rr_n3a3));
            a1a2 = sqrt(dot(rr_a1a2, rr_a1a2));
            a2a3 = sqrt(dot(rr_a2a3, rr_a2a3));

            logger.info(String.format("na: n2a2, n3a3 = %9.3f%9.3f%9.3f%9.3f\n", len0[2], n2a2, len0[5], n3a3));
            logger.info(String.format("ac: a1c1, a2c2 = %9.3f%9.3f%9.3f%9.3f\n", len0[0], a1c1, len0[3], a2c2));
            logger.info(String.format("cn: c1n2, c2n3 = %9.3f%9.3f%9.3f%9.3f\n", len0[1], c1n2, len0[4], c2n3));
            logger.info(String.format("aa: a1a2, a2a3 = %9.3f%9.3f%9.3f%9.3f\n", len_aa[1], a1a2, len_aa[2], a2a3));

            for (i = 0; i < 3; i++) {
                tmp_array[i] = rr_a1a2[i] / a1a2;
            }
            calc_bnd_ang(b_a1a3, tmp_array, a3a1a2);

            for (i = 0; i < 3; i++) {
                tmp_array[i] = rr_a2a3[i] / a2a3;
            }
            calc_bnd_ang(tmp_array, b_a1a3, a2a3a1);

            for (i = 0; i < 3; i++) {
                tmp_array[i] = rr_a1c1[i] / a1c1;
            }
            calc_bnd_ang(b_a1n1, tmp_array, n1a1c1);

            for (i = 0; i < 3; i++) {
                tmp_array[i] = -rr_n3a3[i] / n3a3;
            }
            calc_bnd_ang(b_a3c3, tmp_array, n3a3c3);

            for (i = 0; i < 3; i++) {
                tmp_array1[i] = rr_a2c2[i] / a2c2;
                tmp_array2[i] = -rr_n2a2[i] / n2a2;
            }

            calc_bnd_ang(tmp_array1, tmp_array2, n2a2c2);

            logger.info(String.format("ang_nac = %9.3f%9.3f\n", toDegrees(b_ang0[0]), toDegrees(n1a1c1[0])));
            logger.info(String.format("ang_nac = %9.3f%9.3f\n", toDegrees(b_ang0[3]), toDegrees(n2a2c2[0])));
            logger.info(String.format("ang_nac = %9.3f%9.3f\n", toDegrees(b_ang0[6]), toDegrees(n3a3c3[0])));

            for (i = 0; i < 3; i++) {
                tmp_array1[i] = rr_a1c1[i] / a1c1;
                tmp_array2[i] = rr_c1n2[i] / c1n2;
                tmp_array3[i] = rr_n2a2[i] / n2a2;
            }

            calc_dih_ang(tmp_array1, tmp_array2, tmp_array3, a1c1n2a2);

            for (i = 0; i < 3; i++) {
                tmp_array1[i] = rr_a2c2[i] / a2c2;
                tmp_array2[i] = rr_c2n3[i] / c2n3;
                tmp_array3[i] = rr_n3a3[i] / n3a3;
            }

            calc_dih_ang(tmp_array1, tmp_array2, tmp_array3, a2c2n3a3);

            logger.info(String.format("t_ang1 = %9.3f%9.3f\n", toDegrees(t_ang0[0]), toDegrees(a1c1n2a2[0])));
            logger.info(String.format("t_ang2 = %9.3f%9.3f\n", toDegrees(t_ang0[1]), toDegrees(a2c2n3a3[0])));
        }
    }

    public double calc_t2(double t0) {
        double tmp_value;
        double B0, B1, B2, A0, A1, A2, A3, A4, B2_2, B2_3;
        double K0, K1, K2, K3, t0_2, t0_3, t0_4;

        t0_2 = t0 * t0;
        t0_3 = t0_2 * t0;
        t0_4 = t0_3 * t0;

        A0 = Q[0][0] + Q[0][1] * t0 + Q[0][2] * t0_2 + Q[0][3] * t0_3 + Q[0][4] * t0_4;
        A1 = Q[1][0] + Q[1][1] * t0 + Q[1][2] * t0_2 + Q[1][3] * t0_3 + Q[1][4] * t0_4;
        A2 = Q[2][0] + Q[2][1] * t0 + Q[2][2] * t0_2 + Q[2][3] * t0_3 + Q[2][4] * t0_4;
        A3 = Q[3][0] + Q[3][1] * t0 + Q[3][2] * t0_2 + Q[3][3] * t0_3 + Q[3][4] * t0_4;
        A4 = Q[4][0] + Q[4][1] * t0 + Q[4][2] * t0_2 + Q[4][3] * t0_3 + Q[4][4] * t0_4;

        B0 = R[0][0] + R[0][1] * t0 + R[0][2] * t0_2;
        B1 = R[1][0] + R[1][1] * t0 + R[1][2] * t0_2;
        B2 = R[2][0] + R[2][1] * t0 + R[2][2] * t0_2;

        B2_2 = B2 * B2;
        B2_3 = B2_2 * B2;

        K0 = A2 * B2 - A4 * B0;
        K1 = A3 * B2 - A4 * B1;
        K2 = A1 * B2_2 - K1 * B0;
        K3 = K0 * B2 - K1 * B1;
        tmp_value = (K3 * B0 - A0 * B2_3) / (K2 * B2 - K3 * B1);

        return tmp_value;
    }

    public double calc_t1(double t0, double t2) {
        double tmp_value;
        double U11, U12, U13, U31, U32, U33;
        double t0_2, t2_2;

        t0_2 = t0 * t0;
        t2_2 = t2 * t2;
        U11 = C0[0][0] + C0[0][1] * t0 + C0[0][2] * t0_2;
        U12 = C1[0][0] + C1[0][1] * t0 + C1[0][2] * t0_2;
        U13 = C2[0][0] + C2[0][1] * t0 + C2[0][2] * t0_2;
        U31 = C0[1][0] + C0[1][1] * t2 + C0[1][2] * t2_2;
        U32 = C1[1][0] + C1[1][1] * t2 + C1[1][2] * t2_2;
        U33 = C2[1][0] + C2[1][1] * t2 + C2[1][2] * t2_2;

        tmp_value = (U31 * U13 - U11 * U33) / (U12 * U33 - U13 * U32);

        return (tmp_value);

    }

    public void calc_dih_ang(double[] r1, double[] r2, double[] r3, double[] angle) {
        double[] p = new double[3];
        double[] q = new double[3];
        double[] s = new double[3];
        double arg;

        VectorMath.cross(r1, r2, p);
        VectorMath.cross(r2, r3, q);
        VectorMath.cross(r3, r1, s);

        arg = dot(p, q) / sqrt(dot(p, p) * dot(q, q));
        arg = sign(FastMath.min(FastMath.abs(arg), 1.0e0), arg);
        angle[0] = sign(acos(arg), dot(s, r2));
        //System.out.println("angle_dih: "+angle[0]+"\n");

    }

    public void calc_bnd_ang(double[] r1, double[] r2, double[] angle) {
        double arg;

        arg = dot(r1, r2);
        arg = sign(FastMath.min(FastMath.abs(arg), 1.0e0), arg);
        angle[0] = acos(arg);
    }

    //quotient of two vectors in three dimensional space
    public void quaternion(double[] axis, double quarter_ang, double[] p) {
        double tan_w, tan_sqr, tan1, cosine, sine;

        tan_w = FastMath.tan(quarter_ang);
        tan_sqr = tan_w * tan_w;
        tan1 = 1.0e0 + tan_sqr;
        cosine = (1.0e0 - tan_sqr) / tan1;
        sine = 2.0e0 * tan_w / tan1;

        p[0] = cosine;
        p[1] = axis[0] * sine;
        p[2] = axis[1] * sine;
        p[3] = axis[2] * sine;

        //print statement for debugging purposes
        //System.out.println("Tan_w: " +tan_w +"\n"+"Tan_sqr: " +tan_sqr +"\n"
        // +"Tan1: " +tan1 +"\n"+ "cosine: " +cosine +"\n"+"Sine: " +sine +"\n"
        // +"P: " +p[0] +"\n"+"p: " +p[1] +"\n"+"p: " +p[2] +"\n"+"p: " +p[3] 
        // +"\n");

    }

    //means of representing a rotation of an axis about an origin
    //      with no angular specification
    public void rotation_matrix(double[] q, double[][] U) 
    {
        double q0, q1, q2, q3, b0, b1, b2, b3, q00, q01, q02, q03, q11, q12, q13, q22, q23, q33;

        q0 = q[0];
        q1 = q[1];
        q2 = q[2];
        q3 = q[3];
        //System.out.println("q0: "+q0+" "+q1+" "+q2+" "+q3+"\n");

        b0 = 2.0e0 * q0;
        b1 = 2.0e0 * q1;
        //System.out.println("b0: "+b0+" "+b1+"\n");


        q00 = b0 * q0 - 1.0e0;
        q02 = b0 * q2;
        q03 = b0 * q3;
        //System.out.println("q00: "+q00+" "+q02+" "+q03+" "+"\n");

        q11 = b1 * q1;
        q12 = b1 * q2;
        q13 = b1 * q3;
        //System.out.println("q11: "+q11+" "+q12+" "+q13+"\n");

        b2 = 2.0e0 * q2;
        b3 = 2.0e0 * q3;
        //System.out.println("b2: "+b2+" "+b3+"\n");

        q01 = b0 * q1;
        q22 = b2 * q2;
        q23 = b2 * q3;
        q33 = b3 * q3;
        //System.out.println("q01: "+q01+" "+q22+" "+q23+" "+q33+"\n");

        U[0][0] = q00 + q11;
        U[0][1] = q12 - q03;
        U[0][2] = q13 + q02;
        //System.out.println("U0: "+U[0][0]+" "+U[0][1]+" "+U[0][2]+"\n");

        U[1][0] = q12 + q03;
        U[1][1] = q00 + q22;
        U[1][2] = q23 - q01;
        //System.out.println("U1: "+U[1][0]+" "+U[1][1]+" "+U[1][2]+"\n");


        U[2][0] = q13 - q02;
        U[2][1] = q23 + q01;
        U[2][2] = q00 + q33;
        //System.out.println("U2: "+U[2][0]+" "+U[2][1]+" "+U[2][2]+"\n");
    }
}
