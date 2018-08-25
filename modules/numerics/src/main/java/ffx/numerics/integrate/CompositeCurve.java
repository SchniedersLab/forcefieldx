/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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

import java.util.Arrays;
import java.util.List;

/**
 * A CompositeCurve represents points along a sum of functions which also extend
 * FunctionDataCurve.
 *
 * @author Jacob M. Litman
 */
public class CompositeCurve extends FunctionDataCurve {
    private final FunctionDataCurve[] curves;
    private final double[] coeffs;
    private final int nCurves;

    /**
     * Constructs a CompositeCurve that aggregates multiple FunctionDataCurves
     * with variable weights to each component FunctionDataCurve.
     *
     * @param componentCurves Underlying FunctionDataCurves
     * @param coefficients Weight to each component curve
     */
    public CompositeCurve(List<FunctionDataCurve> componentCurves, List<Double> coefficients) {
        assert !componentCurves.isEmpty();
        nCurves = componentCurves.size();
        this.curves = new FunctionDataCurve[nCurves];
        componentCurves.toArray(this.curves);
        assert (coefficients == null || nCurves == coefficients.size());
        
        if (coefficients == null) {
            coeffs = new double[nCurves];
            Arrays.fill(coeffs, 1.0);
        } else {
            coeffs = coefficients.stream().mapToDouble(Double::doubleValue).toArray();
        }
        
        FunctionDataCurve curve0 = curves[0];
        lb = curve0.lowerBound();
        ub = curve0.upperBound();
        halfWidthEnd = curve0.halfWidthEnds();
        x = curve0.getX();
        double sep = curve0.binWidth();
        
        boolean isInvalid = Arrays.stream(curves).anyMatch((FunctionDataCurve c) -> {
            if (lb != c.lowerBound()) {
                return true;
            }
            if (ub != c.upperBound()) {
                return true;
            }
            if (halfWidthEnd != c.halfWidthEnds()) {
                return true;
            }
            return sep != c.binWidth();
        });
        if (isInvalid) {
            throw new IllegalArgumentException(" Not all curves passed to CompositeCurve had the same x[] points!");
        }
        
        int nPoints = curve0.numPoints();
        points = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            points[i] = 0.0;
            for (int j = 0; j < nCurves; j++) {
                points[i] += (coeffs[j] * curves[j].getFxPoint(i));
            }
        }
    }

    /** {@inheritDoc} */
    @Override
    public double integralAt(double x) {
        double val = 0.0;
        for (int i = 0; i < nCurves; i++) {
            val += (curves[i].integralAt(x) * coeffs[i]);
        }
        return val;
    }
    
    /** {@inheritDoc} */
    @Override
    public double fX(double x) {
        return valAt(x);
    }
    
    // Private, non-overrideable method for use in the constructor.
    private double valAt(double x) {
        double val = 0.0;
        for (int i = 0; i < nCurves; i++) {
            val += (curves[i].fX(x) * coeffs[i]);
        }
        return val;
    }

    /**
     * Gets the component FunctionDataCurves of this CompositeCurve.
     *
     * @return List of component FunctionDataCurves.
     */
    public List<FunctionDataCurve> getSubCurves() {
        return Arrays.asList(curves);
    }

    /**
     * Gets the weights to the corresponding component curves.
     *
     * @return Constant weights
     */
    public double[] getWeights() {
        return Arrays.copyOf(coeffs, coeffs.length);
    }
    
    /** {@inheritDoc} */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(String.format("Composite curve with %d points from lower bound %9.3g and upper bound %9.3g", points.length, lb, ub));
        if (halfWidthEnd) {
            sb.append(" and half-width start/end bins.\nComponent curves:\n");
        }
        for (FunctionDataCurve curve : curves) {
            sb.append(curve.toString());
        }
        return sb.toString();
    }
}
