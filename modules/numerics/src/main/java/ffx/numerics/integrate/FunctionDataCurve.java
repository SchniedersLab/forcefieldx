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
package ffx.numerics.integrate;

/**
 * A FunctionDataCurve represents a set of points along a 1-dimensional, 
 * analytically integrable function.
 * @author Jacob M. Litman
 */
public abstract class FunctionDataCurve implements DataSet {
    protected double lb;
    protected double ub;
    protected double[] points;
    protected boolean halfWidthEnd;
    protected double[] x;
    
    @Override
    public double lowerBound() {
        return lb;
    }

    @Override
    public double upperBound() {
        return ub;
    }
    
    @Override
    public int numPoints() {
        return points.length;
    }
    
    @Override
    public double binWidth() {
        double divisor = halfWidthEnds() ? (double) (points.length - 2) : (double) (points.length - 1);
        return (ub - lb) / divisor;
    }
    
    /**
     * Evaluates the functions analytical integral over the entire range of points.
     * @return Exact finite integral
     */
    public double analyticalIntegral() {
        return anaylticalIntegral(lowerBound(), upperBound());
    }
    
    /**
     * Evaluates the function's analytical integral over a range.
     * @param lb Lower integration bound
     * @param ub Upper integration bound
     * @return Exact finite integral of range
     */
    public double anaylticalIntegral(double lb, double ub) {
        return integralAt(ub) - integralAt(lb);
    }
    
    /**
     * Analytical integral at a point.
     * @param x Point
     * @return Exact finite integral of 0 to this point
     */
    public abstract double integralAt(double x);
    
    /**
     * Evaluates the function at x.
     * @param x x
     * @return f(x)
     */
    public abstract double fX(double x);
    
    @Override
    public double getFxPoint(int index) {
        return points[index];
    }
    
    @Override
    public double[] getAllFxPoints() {
        int npoints = points.length;
        double[] retArray = new double[npoints];
        System.arraycopy(points, 0, retArray, 0, npoints);
        return retArray;
    }
    
    /**
     * Used to check that the passed-in x array is composed of equally-spaced
     * points from lb to ub.
     * @param x 
     */
    protected final void assertXIntegrity(double[] x) {
        assert ub > lb;
        int nX = numPoints();
        double sep = binWidth();
        if (halfWidthEnd) {
            assert x.length == nX;
            assert lb == x[0];
            assert ub == x[nX-1];
            
            assert approxEquals(x[1], lb + 0.5*sep);
            assert approxEquals(x[nX-2], (ub - 0.5*sep));
            
            for (int i = 2; i < (nX - 2); i++) {
                double target = lb + 0.5 * sep;
                target += ((i-1) * sep);
                assert approxEquals(x[i], target);
            }
        } else {
            for (int i = 0; i < x.length; i++) {
                assert approxEquals(x[i], x[0] + i * sep);
            }
        }
    }
    
    /**
     * Checks for equality to +/- 10 ulp.
     * @param x1
     * @param x2
     * @return 
     */
    public static boolean approxEquals(double x1, double x2) {
        return (approxEquals(x1, x2, 10.0));
    }
    
    /**
     * Compare two doubles to machine precision.
     * @param x1
     * @param x2
     * @param ulpMult
     * @return 
     */
    public static boolean approxEquals(double x1, double x2, double ulpMult) {
        double diff = Math.abs(x1 - x2);
        // Ulp the larger of the two values.
        double ulp = Math.ulp(Math.max(Math.abs(x1), Math.abs(x2)));
        return (diff < (ulp * ulpMult));
    }
    
    @Override
    public double[] getX() {
        double[] copyX = new double[x.length];
        System.arraycopy(x, 0, copyX, 0, x.length);
        return copyX;
    }
    
    @Override
    public boolean halfWidthEnds() {
        return halfWidthEnd;
    }
    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(String.format("Function f(x) curve with %d points from lower bound %9.3g and upper bound %9.3g", points.length, lb, ub));
        if (halfWidthEnd) {
            sb.append(" and half-width start/end bins");
        }
        sb.append(".");
        return sb.toString();
    }
}
