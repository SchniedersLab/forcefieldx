// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// ******************************************************************************
package ffx.potential.utils;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.asin;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.atan;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.signum;


import ffx.potential.bonded.Atom;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * Inertia computes the principal moments of inertia for the system, and optionally translates the
 * center of mass to the origin and rotates the principal axes onto the global axes.
 *
 * Reference:
 * Herbert Goldstein, "Classical Mechanics, 2nd Edition",
 *  Addison-Wesley, Reading, MA, 1980; see the Euler angle
 *  xyz convention in Appendix B
 *
 * @author Aaron J. Nessler
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Inertia {

    private static final Logger logger = Logger.getLogger(Inertia.class.getName());

    /**
     * Compute the moments of inertia for all atoms in the supplied array.
     *
     * @param atoms Atom array.
     * @return The moments of inertia.
     */
    public static double[][] momentsOfInertia(Atom[] atoms, boolean moved, boolean print) {
        double[] mass = new double[atoms.length];
        int nAtoms = atoms.length;
        double[] x = new double[nAtoms];
        double[] y = new double[nAtoms];
        double[] z = new double[nAtoms];

        int index = 0;
        for (Atom atom : atoms) {
            mass[index] = atom.getMass();
            x[index] = atom.getX();
            y[index] = atom.getY();
            z[index] = atom.getZ();
            index++;
        }

        return momentsOfInertia(x, y, z, mass, moved, print);
    }

    /**
     * Compute the moments of inertia
     *
     * @param x Array of atomic X-coordinates.
     * @param y Array of atomic X-coordinates.
     * @param z Array of atomic X-coordinates.
     * @return The moment of inertia.
     */
    public static double[][] momentsOfInertia(double[] x, double[] y, double[] z, double[] mass, boolean moved, boolean print) {
        assert (x.length == y.length);
        assert (y.length == z.length);

        // Find the centroid of the atomic coordinates.
        double total = 0.0;
        double xcm = 0.0;
        double ycm = 0.0;
        double zcm = 0.0;
        int nAtoms = x.length;
        for (int i = 0; i < nAtoms; i++) {
            double massValue = mass[i];
            total += massValue;
            xcm += x[i] * massValue;
            ycm += y[i] * massValue;
            zcm += z[i] * massValue;
        }
        xcm /= total;
        ycm /= total;
        zcm /= total;

        // Compute and then diagonalize the inertia tensor
        double xx = 0.0;
        double xy = 0.0;
        double xz = 0.0;
        double yy = 0.0;
        double yz = 0.0;
        double zz = 0.0;

        double xterm;
        double yterm;
        double zterm;

        for(int i = 0; i < nAtoms; i++){
            double massValue = mass[i];
            xterm = x[i] - xcm;
            yterm = y[i] - ycm;
            zterm = z[i] - zcm;
            xx += xterm * xterm * massValue;
            xy += xterm * yterm * massValue;
            xz += xterm * zterm * massValue;
            yy += yterm * yterm * massValue;
            yz += yterm * zterm * massValue;
            zz += zterm * zterm * massValue;
        }
        double[][] tensor = new double[3][3];
        tensor[0][0] = yy + zz;
        tensor[0][1] = -xy;
        tensor[0][2] = -xz;
        tensor[1][0] = -xy;
        tensor[1][1] = xx + zz;
        tensor[1][2] = -yz;
        tensor[2][0] = -xz;
        tensor[2][1] = -yz;
        tensor[2][2] = xx + yy;

        // Diagonalize the matrix
        Array2DRowRealMatrix cMatrix = new Array2DRowRealMatrix(tensor, false);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(cMatrix);
        // Extract the quaternions.
        double[] moment = eigenDecomposition.getRealEigenvalues();
        Arrays.sort(moment);
        // Initialization below to match Tinker inertia.f
        double[][] vec = new double[3][3];
        vec[2] = eigenDecomposition.getEigenvector(0).toArray();
        vec[1] = eigenDecomposition.getEigenvector(1).toArray();
        vec[0] = eigenDecomposition.getEigenvector(2).toArray();
//        double[][] vec = eigenDecomposition.getVT().getData();

        // Select the direction for each principal moment axis
        double dot = 0.0;
        for(int i = 0; i < 2; i++){
            for(int j = 0; j < nAtoms; j++){
                xterm = vec[i][0] * (x[j] - xcm);
                yterm = vec[i][1] * (y[j] - ycm);
                zterm = vec[i][2] * (z[j] - zcm);
                dot = xterm + yterm + zterm;
                if(dot < 0.0){
                    for(int k = 0; k < 3; k++){
                        vec[i][k] = -vec[i][k];
                    }
                }
                if(dot != 0.0){
                    break;
                }
            }
            if(dot != 0.0){
                break;
            }
        }

        // moment axes must give a right-handed coordinate system
        xterm = vec[0][0] * (vec[1][1] * vec[2][2] - vec[2][1] * vec[1][2]);
        yterm = vec[0][1] * (vec[2][0] * vec[1][2] - vec[1][0] * vec[2][2]);
        zterm = vec[0][2] * (vec[1][0] * vec[2][1] - vec[2][0] * vec[1][1]);
        dot = xterm + yterm + zterm;
        if(dot < 0.0){
            for(int i = 0; i < 3; i++){
                vec[2][i] = -vec[2][i];
            }
        }

        // principal moment axes form rows of Euler rotation matrix
        if(moved){
            double[][] a = new double[3][3];
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; j++){
                    a[j][i] = vec[j][i];
                }
            }
            // translate to origin, then apply Euler rotation matrix
            for(int i = 0; i < nAtoms; i++){
                xterm = x[i] - xcm;
                yterm = y[i] - ycm;
                zterm = z[i] - zcm;
                x[i] = a[0][0] * xterm + a[1][0] * yterm + a[2][0] * zterm;
                y[i] = a[0][1] * xterm + a[1][1] * yterm + a[2][1] * zterm;
                z[i] = a[0][2] * xterm + a[1][2] * yterm + a[2][2] * zterm;
            }
        }

        // print the center of mass and Euler angle values
        if(print){
            logger.info(format("\n Center of Mass Coordinates: %8.4f %8.4f %8.4f", xcm, ycm, zcm));
            // invert vec
            RealMatrix m = new Array2DRowRealMatrix(vec, true);
            m = new LUDecomposition(m).getSolver().getInverse();
            vec = m.getData();
            double radian = 180/PI;
            double[] angles = roteuler(vec);
            // Convert to degrees
            for(int i = 0; i< 3; i++){
                angles[i] *= radian;
            }
            logger.info(format(" Euler Angles (Phi/Theta/Psi): %8.3f %8.3f %8.3f", angles[0], angles[1], angles[2]));
            logger.info(" Moments of Inertia and Principle Axes:\n  Moments (amu Ang^2): \t X-, Y-, and Z-Components of Axes:");
            for(int i = 0; i < 3; i++){
                logger.info(format("  %16.3f %12.6f %12.6f %12.6f", moment[i], vec[i][0], vec[i][1], vec[i][2]));
            }
        }

        double[][] momentsAndVectors = new double[3][4];
        for(int i = 0; i< 3;i++){
            for(int j = 0; j < 4; j++){
                if(j==0) {
                    momentsAndVectors[i][j] = moment[i];
                }else{
                    momentsAndVectors[i][j] = vec[i][j-1];
                }
            }
        }

        return momentsAndVectors;
    }

    /**
     * Compute the moments of inertia.
     *
     * @param xyz Array of atomic coordinates (xyz = [X0, Y0, Z0, X1, Y1, Z1, ...].
     * @return The radius of gyration.
     */
    public static double[][] momentsOfInertia(double[] xyz, double[] mass, boolean moved, boolean print) {
        assert (xyz.length % 3 == 0);
        int nAtoms = xyz.length / 3;
        // Find the centroid of the atomic coordinates.
        double[] x = new double[nAtoms];
        double[] y = new double[nAtoms];
        double[] z = new double[nAtoms];

        for (int i = 0; i < nAtoms; i++) {
            int index = i * 3;
            x[i] = xyz[index++];
            y[i] = xyz[index++];
            z[i] = xyz[index];
        }

        return momentsOfInertia(x, y, z, mass, moved, print);
    }

    /**
     * Computes a set of Euler angle values consistent with an input rotation matrix
     */
    public static double[] roteuler(double[][] a){
        // set the tolerance for Euler angles and rotation elements
        double eps = 1.0E-7;

        //  get a trial value of theta from a single rotation element
        double psi = 0.0;
        double phi = 0.0;
        double theta = asin(min(1.0, max(-1.0, -a[2][0])));
        double ctheta = cos(theta);
        double stheta = -a[2][0];

        if(abs(ctheta) <= eps){
            // set the phi/psi difference when theta is either 90 or -90
            if(abs(a[0][2]) < eps){
                psi = asin(min(1.0, max(-1.0, -a[0][1]/a[2][0])));
            }else if(abs(a[0][1]) < eps){
                psi = acos(min(1.0, max(-1.0, -a[0][2]/a[2][0])));
            }else{
                psi = atan(a[0][1]/a[0][2]);
            }
        }else{
            // set the phi and psi values for all other theta values
            if(abs(a[0][0]) < eps){
                phi = asin(min(1.0,max(-1.0,a[1][0]/ctheta)));
            }else if(abs(a[1][0]) < eps){
                phi = acos(min(1.0,max(-1.0,a[0][0]/ctheta)));
            }else{
                phi = atan(a[1][0]/a[0][0]);
            }
            if(abs(a[2][2]) < eps){
                psi = asin(min(1.0,max(-1.0, a[2][1]/ctheta)));
            }else if(abs(a[2][1]) < eps){
                psi = acos(min(1.0,max(-1.0,a[2][2]/ctheta)));
            }else{
                psi = atan(a[2][1] / a[2][2]);
            }
        }

        // find sine and cosine of the trial phi and psi values
        double cphi = cos(phi);
        double sphi = sin(phi);
        double cpsi = cos(psi);
        double spsi = sin(psi);

        //  reconstruct the diagonal of the rotation matrix'
        double[] b = new double[3];
        b[0] = ctheta * cphi;
        b[1] = spsi * stheta * sphi + cpsi * cphi;
        b[2] = ctheta * cpsi;

        // compare the correct matrix diagonal to rebuilt diagonal
        boolean[] flip = new boolean[3];
        for(int i = 0; i < 3; i++){
            flip[i] = abs(a[i][i]-b[i]) > eps;
        }
        // alter Euler angles to get correct rotation matrix values

        phi = (flip[0] && flip[1]) ? phi - PI * signum(phi): phi;
        theta = (flip[0] && flip[2]) ? -theta + PI * signum(theta): theta;
        psi = (flip[1] && flip[2]) ? psi - PI * signum(psi): psi;

        // convert maximum negative angles to positive values
        phi = (phi <= -PI) ? PI:phi;
        theta = (theta <= -PI) ? PI:theta;
        psi = (psi <= -PI) ? PI:psi;

        return new double[]{phi, theta, psi};
    }
}

