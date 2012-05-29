/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012
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
package ffx.autoparm;

import java.io.*;
import java.text.DecimalFormat;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

/**
 * Superpose two molecules.
 *
 * @author Gaurav Chattree
 * @since 1.0
 * @version $Id: $
 */
public class Superpose {

    int n1, n2, nfit;
    double xyz1[][], xyz2[][];
    String name1[], name2[];
    double xmid = 0, ymid = 0, zmid = 0;

    /**
     * <p>Constructor for Superpose.</p>
     *
     * @param xyzfname1 a {@link java.lang.String} object.
     * @param xyzfname2 a {@link java.lang.String} object.
     */
    public Superpose(String xyzfname1, String xyzfname2) {

        File structure_xyz1 = new File(xyzfname1);
        int n = 1;
        String oxyz1fname = null;
        String old1 = xyzfname1;
        while (structure_xyz1 != null && structure_xyz1.exists() && structure_xyz1.canRead()) {
            oxyz1fname = xyzfname1;
            n++;
            xyzfname1 = old1;
            xyzfname1 = xyzfname1 + "_" + Integer.toString(n);
            structure_xyz1 = new File(xyzfname1);
        }
        structure_xyz1 = new File(oxyz1fname);

        File structure_xyz2 = new File(xyzfname2);
        n = 1;
        String oxyz2fname = null;
        String old2 = xyzfname2;
        while (structure_xyz2 != null && structure_xyz2.exists() && structure_xyz2.canRead()) {
            oxyz2fname = xyzfname2;
            n++;
            xyzfname2 = old2;
            xyzfname2 = xyzfname2 + "_" + Integer.toString(n);
            structure_xyz2 = new File(xyzfname2);
        }
        structure_xyz2 = new File(oxyz2fname);

        String line;
        String[] ln;

        //Create two arrays with the xyz coordinates for each atom in each molecule
        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(structure_xyz1)));
            line = br.readLine();
            n1 = Integer.parseInt(line.trim());
            name1 = new String[n1];
            xyz1 = new double[n1][3];
            int i = 0;
            while ((line = br.readLine()) != null) {
                ln = line.trim().split(" +");
                name1[i] = ln[1];
                xyz1[i][0] = Double.parseDouble(ln[2]);
                xyz1[i][1] = Double.parseDouble(ln[3]);
                xyz1[i][2] = Double.parseDouble(ln[4]);
                i++;
            }
            br.close();
            br = new BufferedReader(new InputStreamReader(new FileInputStream(structure_xyz2)));
            line = br.readLine();
            n2 = Integer.parseInt(line.trim());
            name2 = new String[n2];
            xyz2 = new double[n2][3];
            i = 0;
            while ((line = br.readLine()) != null) {
                ln = line.trim().split(" +");
                name2[i] = ln[1];
                xyz2[i][0] = Double.parseDouble(ln[2]);
                xyz2[i][1] = Double.parseDouble(ln[3]);
                xyz2[i][2] = Double.parseDouble(ln[4]);
                i++;
            }
            br.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        nfit = Math.min(n1, n2);


        center();
        quatfit();

        //The first set of coordinates is returned to its original position.
        for (int i = 0; i < n1; i++) {
            xyz1[i][0] = xyz1[i][0] + xmid;
            xyz1[i][1] = xyz1[i][1] + ymid;
            xyz1[i][2] = xyz1[i][2] + zmid;
        }
        for (int i = 0; i < n2; i++) {
            xyz2[i][0] = xyz2[i][0] + xmid;
            xyz2[i][1] = xyz2[i][1] + ymid;
            xyz2[i][2] = xyz2[i][2] + zmid;
        }
        output();
    }

    /**
     * <p>center</p>
     */
    public void center() {
        double norm = 0;
        for (int i = 0; i < n2; i++) {
            xmid += xyz2[i][0];
            ymid += xyz2[i][1];
            zmid += xyz2[i][2];
            norm++;
        }
        xmid = xmid / norm;
        ymid = ymid / norm;
        zmid = zmid / norm;
        for (int i = 0; i < n2; i++) {
            xyz2[i][0] = xyz2[i][0] - xmid;
            xyz2[i][1] = xyz2[i][1] - ymid;
            xyz2[i][2] = xyz2[i][2] - zmid;
        }
        xmid = 0;
        ymid = 0;
        zmid = 0;
        norm = 0;
        for (int i = 0; i < n1; i++) {
            xmid += xyz1[i][0];
            ymid += xyz1[i][1];
            zmid += xyz1[i][2];
            norm++;
        }
        xmid = xmid / norm;
        ymid = ymid / norm;
        zmid = zmid / norm;
        for (int i = 0; i < n1; i++) {
            xyz1[i][0] = xyz1[i][0] - xmid;
            xyz1[i][1] = xyz1[i][1] - ymid;
            xyz1[i][2] = xyz1[i][2] - zmid;
        }
    }

    /**
     * <p>quatfit</p>
     */
    public void quatfit() {
        double xxyx = 0.0;
        double xxyy = 0.0;
        double xxyz = 0.0;
        double xyyx = 0.0;
        double xyyy = 0.0;
        double xyyz = 0.0;
        double xzyx = 0.0;
        double xzyy = 0.0;
        double xzyz = 0.0;
        double c[][] = new double[4][4];
        for (int i = 0; i < nfit; i++) {
            xxyx = xxyx + xyz1[i][0] * xyz2[i][0];
            xxyy = xxyy + xyz1[i][1] * xyz2[i][0];
            xxyz = xxyz + xyz1[i][2] * xyz2[i][0];
            xyyx = xyyx + xyz1[i][0] * xyz2[i][1];
            xyyy = xyyy + xyz1[i][1] * xyz2[i][1];
            xyyz = xyyz + xyz1[i][2] * xyz2[i][1];
            xzyx = xzyx + xyz1[i][0] * xyz2[i][2];
            xzyy = xzyy + xyz1[i][1] * xyz2[i][2];
            xzyz = xzyz + xyz1[i][2] * xyz2[i][2];
        }
        c[0][0] = xxyx + xyyy + xzyz;
        c[0][1] = xzyy - xyyz;
        c[1][1] = xxyx - xyyy - xzyz;
        c[0][2] = xxyz - xzyx;
        c[1][2] = xxyy + xyyx;
        c[2][2] = xyyy - xzyz - xxyx;
        c[0][3] = xyyx - xxyy;
        c[1][3] = xzyx + xxyz;
        c[2][3] = xyyz + xzyy;
        c[3][3] = xzyz - xxyx - xyyy;

        RealMatrix a = new Array2DRowRealMatrix(new double[][]{
                    {c[0][0], c[0][1], c[0][2], c[0][3]},
                    {c[0][1], c[1][1], c[1][2], c[1][3]}, {c[0][2], c[1][2], c[2][2], c[2][3]},
                    {c[0][3], c[1][3], c[2][3], c[3][3]}});
        EigenDecompositionImpl e = new EigenDecompositionImpl(a, 1);
        a = e.getV();
        double[] q = {a.getEntry(0, 0), a.getEntry(1, 0), a.getEntry(2, 0), a.getEntry(3, 0)};
        double rot[][] = new double[3][3];
        rot[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
        rot[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
        rot[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);
        rot[0][1] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
        rot[1][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
        rot[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);
        rot[0][2] = 2.0 * (q[3] * q[1] - q[0] * q[2]);
        rot[1][2] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
        rot[2][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
        double xrot, yrot, zrot;
        for (int i = 0; i < n2; i++) {
            xrot = xyz2[i][0] * rot[0][0] + xyz2[i][1] * rot[0][1] + xyz2[i][2] * rot[0][2];
            yrot = xyz2[i][0] * rot[1][0] + xyz2[i][1] * rot[1][1] + xyz2[i][2] * rot[1][2];
            zrot = xyz2[i][0] * rot[2][0] + xyz2[i][1] * rot[2][1] + xyz2[i][2] * rot[2][2];
            xyz2[i][0] = xrot;
            xyz2[i][1] = yrot;
            xyz2[i][2] = zrot;
        }
    }

    /**
     * <p>rms</p>
     *
     * @return a double.
     */
    public double rms() {
        double xr = 0, yr = 0, zr = 0, rms = 0, norm = 0;
        for (int i = 0; i < nfit; i++) {
            xr = xyz1[i][0] - xyz2[i][0];
            yr = xyz1[i][1] - xyz2[i][1];
            zr = xyz1[i][2] - xyz2[i][2];
            rms += xr * xr + yr * yr + zr * zr;
            norm++;
        }
        rms = Math.sqrt(rms / norm);
        return rms;
    }

    /**
     * <p>output</p>
     */
    public void output() {
        DecimalFormat myFormatter = new DecimalFormat(" ##########0.000000;-##########0.000000");
        String headings = String.format("\n   %9s%25s%20s\n  %9s%26s%17s\n", "Atom in the", "Atom in the", "Distance", "First Structure", "Second Structure", "Separated");
        System.out.println(headings);
        double xr, yr, zr, dist;
        String out;
        for (int i = 0; i < nfit; i++) {
            xr = xyz1[i][0] - xyz2[i][0];
            yr = xyz1[i][1] - xyz2[i][1];
            zr = xyz1[i][2] - xyz2[i][2];
            dist = Math.sqrt(xr * xr + yr * yr + zr * zr);
            out = String.format("%8d-%s%24d-%s%23s", (i + 1), name1[i], (i + 1), name2[i], myFormatter.format(dist));
            System.out.println(out);
        }
        System.out.println("\n Root Mean Square Distance:                       " + myFormatter.format(rms()));

    }

    /**
     * <p>main</p>
     *
     * @param args an array of {@link java.lang.String} objects.
     */
    public static void main(String args[]) {
        Superpose s = new Superpose("/home/gchattree/Research/Compounds/s_test3_compounds/famotidine/famotidine.xyz", "/home/gchattree/Research/Compounds/s_test3_compounds/famotidine/ttt.xyz");
    }
}
