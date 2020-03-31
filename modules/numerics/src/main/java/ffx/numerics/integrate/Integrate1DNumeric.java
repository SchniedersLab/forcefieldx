//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
//******************************************************************************
package ffx.numerics.integrate;

import java.util.logging.Logger;
import java.util.stream.IntStream;

/**
 * This program integrates using four methods: rectangular integration, the
 * trapezoidal method, Simpson's Three Point Integration, and Boole's Five Point
 * Integration. Of these, Simpson's and Boole's integration methods generally
 * perform the best. Parallelized methods exist and seem to function, but the
 * performance of a parallelized method should not be necessary, and the
 * parallelization introduces some complexity and inefficiency.
 * <p>
 * Primary methods: integrateData and integrateByBins, using some DataSet,
 * which is often a DoublesDataSet.
 *
 * @author Claire O'Connell
 * @author Jacob M. Litman
 */
public class Integrate1DNumeric {

    private static final Logger logger = Logger.getLogger(Integrate1DNumeric.class.getName());
    /**
     * Constant for used for Simpson's rule.
     */
    private static final double ONE_THIRD = (1.0 / 3.0);
    /**
     * Constant for used for Boole's rule.
     */
    private static final double BOOLE_FACTOR = (2.0 / 45.0);

    /**
     * Numerically integrates a data set using Boole's rule.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @return Area of integral
     */
    public static double booles(DataSet data, IntegrationSide side) {
        double area = 0;
        int lb = 0;
        int ub = data.numPoints() - 1;
        if (data.halfWidthEnds()) {
            area = trapezoidalEnds(data, side);
            lb++;
            ub--;
        }
        area += booles(data, side, lb, ub);
        return area;
    }

    /**
     * Numerically integrates a data set, in bounds lb-ub inclusive, using
     * Boole's rule.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @param lb   First index to integrate over
     * @param ub   Last index to integrate over
     * @return Area of integral
     */
    public static double booles(DataSet data, IntegrationSide side, int lb, int ub) {
        double area = 0;
        double width = data.binWidth();
        double[] points = data.getAllFxPoints();
        int increment = 4;

        int nBins = (ub - lb) / increment;
        int lowerNeglected;
        int upperNeglected;

        switch (side) {
            case RIGHT:
                for (int i = ub; i > (lb + increment - 1); i -= increment) {
                    area += (7 * points[i] + 32 * points[i - 1] + 12 * points[i - 2] + 32 * points[i - 3] + 7 * points[i - 4]);
                }
                lowerNeglected = lb;
                upperNeglected = ub - (increment * nBins);
                break;
            case LEFT:
            default:
                for (int i = lb; i < (ub - increment + 1); i += increment) {
                    area += (7 * points[i] + 32 * points[i + 1] + 12 * points[i + 2] + 32 * points[i + 3] + 7 * points[i + 4]);
                }
                lowerNeglected = lb + (increment * nBins);
                upperNeglected = ub;
                break;
        }
        area *= BOOLE_FACTOR;
        area *= width;

        area += finishIntegration(data, side, lowerNeglected, upperNeglected, IntegrationType.BOOLE);
        return area;
    }

    /**
     * Numerically integrates a data set using Boole's rule. The sequential
     * version is preferred unless necessary.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @return Area of integral
     */
    public static double boolesParallel(DataSet data, IntegrationSide side) {
        double area = 0;
        int lb = 0;
        int ub = data.numPoints() - 1;
        if (data.halfWidthEnds()) {
            area = trapezoidalEnds(data, side);
            lb++;
            ub--;
        }
        area += boolesParallel(data, side, lb, ub);
        return area;
    }

    /**
     * Numerically integrates a data set, in bounds lb-ub inclusive, using
     * Boole's rule. The sequential version is preferred unless necessary.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @param lb   First index to integrate over
     * @param ub   Last index to integrate over
     * @return Area of integral
     */
    public static double boolesParallel(DataSet data, IntegrationSide side, int lb, int ub) {
        double area = 0;
        double width = data.binWidth();
        double[] points = data.getAllFxPoints();
        int increment = 4;

        int nBins = (ub - lb) / increment;
        int lowerNeglected;
        int upperNeglected;

        switch (side) {
            case RIGHT:
                lowerNeglected = lb;
                upperNeglected = ub - (increment * nBins);
                area = IntStream.range(0, nBins).parallel().mapToDouble((int i) -> {
                    int fromPoint = ub - (increment * i);
                    return (7 * points[fromPoint - 4] + 32 * points[fromPoint - 3] + 12 * points[fromPoint - 2] + 32 * points[fromPoint - 1] + 7 * points[fromPoint]);
                }).sum();
                break;
            case LEFT:
            default:
                lowerNeglected = lb + (increment * nBins);
                upperNeglected = ub;
                area = IntStream.range(0, nBins).parallel().mapToDouble((int i) -> {
                    int fromPoint = lb + (increment * i);
                    return (7 * points[fromPoint] + 32 * points[fromPoint + 1] + 12 * points[fromPoint + 2] + 32 * points[fromPoint + 3] + 7 * points[fromPoint + 4]);
                }).sum();
                break;
        }
        area *= BOOLE_FACTOR;
        area *= width;

        area += finishIntegration(data, side, lowerNeglected, upperNeglected, IntegrationType.BOOLE);
        return area;
    }

    /**
     * Generates a set of points along x.
     *
     * @param lb            Beginning value, inclusive
     * @param ub            Ending value, inclusive
     * @param nPoints       Total number of points
     * @param halfWidthEnds If ends should have 1/2 regular separation
     * @return an array of {@link double} objects.
     */
    public static double[] generateXPoints(double lb, double ub, int nPoints, boolean halfWidthEnds) {
        if (lb >= ub) {
            throw new IllegalArgumentException("ub must be greater than lb");
        }
        double[] points = new double[nPoints];
        int sepDivisor = halfWidthEnds ? nPoints - 2 : nPoints - 1;
        double sep = (ub - lb) / ((double) sepDivisor);

        if (halfWidthEnds) {
            points[0] = lb;
            points[nPoints - 1] = ub;
            for (int i = 1; i < nPoints - 1; i++) {
                points[i] = lb + ((double) i) * sep - 0.5 * sep;
            }
        } else {
            for (int i = 0; i < nPoints; i++) {
                points[i] = lb + ((double) i) * sep;
            }
        }
        return points;
    }

    /**
     * Returns the contribution of each bin to the overall integral as an array;
     * will be most accurate at break-points for the integration type. Overall
     * integral is the sum of the array of doubles, plus or minus minor rounding
     * errors.
     * <p>
     * If N is a breakpoint for Boole's rule:
     * N+1 is a trapezoid from N to N+1.
     * N+2 is Simpson's from N to N+2, minus the prior trapezoid (N+1).
     * N+3 is 4-point integration from N to N+3, minus the N to N+2 parabola.
     * N+4 is the full Boole's Rule from N to N+4, minus the N to N+3 4-point
     * integral.
     *
     * @param data    Data to integrate
     * @param side    Side to integrate from
     * @param maxType Maximum rule to be used
     * @return Per-bin contributions to integral.
     */
    public static double[] integrateByBins(DataSet data, IntegrationSide side, IntegrationType maxType) {
        boolean halfWide = data.halfWidthEnds();
        double[] fX = data.getAllFxPoints();
        int numPoints = data.numPoints();
        double width = data.binWidth();
        double[] vals = new double[numPoints];

        int lb = halfWide ? 1 : 0;
        int ub = halfWide ? numPoints - 2 : numPoints - 1;

        int increment = maxType.pointsNeeded() - 1;
        increment = Math.max(1, increment); // Deal w/ rectangle integration.

        switch (side) {
            case RIGHT:
                vals[ub] = width * fX[ub]; // Begin with rectangle.

                /**
                 * For each bin, its contribution to this sub-window's value is
                 * the integral from (start to bin) minus the integral from
                 * (start to prior bin).
                 */
                for (int i = ub - 1; i >= lb; i--) {
                    int fromUB = ub - i - 1;
                    fromUB /= increment;
                    fromUB *= increment;
                    int lastBegin = ub - fromUB;

                    double val = finishIntegration(data, side, i, lastBegin, maxType);
                    val -= finishIntegration(data, side, i + 1, lastBegin, maxType);
                    vals[i] = val;
                }

                vals[ub - 1] -= vals[ub]; // Remove double-counting at start.

                if (halfWide) {
                    switch (maxType) {
                        case RECTANGULAR:
                            vals[1] += 0.5 * width * fX[1];
                            vals[numPoints - 1] = (0.5 * width * fX[numPoints - 1]);
                            break;
                        default:
                            vals[0] = 0.25 * width * fX[0];
                            vals[1] += 0.25 * width * fX[1];
                            vals[numPoints - 2] += (0.25 * width * fX[numPoints - 2]);
                            vals[numPoints - 1] = (0.25 * width * fX[numPoints - 1]);
                            break;
                    }
                }
                break;
            case LEFT:
            default:
                int shift = halfWide ? 1 : 0;
                vals[lb] = width * fX[lb]; // Begin with rectangle.

                /**
                 * For each bin, its contribution to this sub-window's value is
                 * the integral from (start to bin) minus the integral from
                 * (start to prior bin).
                 */
                for (int i = lb + 1; i <= ub; i++) {
                    // Remove remainder via division-and-multiplication.
                    int lastBegin = ((i - 1 - shift) / increment);
                    lastBegin *= increment;
                    lastBegin += shift;

                    double val = finishIntegration(data, side, lastBegin, i, maxType);
                    val -= finishIntegration(data, side, lastBegin, i - 1, maxType);
                    vals[i] = val;
                }

                vals[lb + 1] -= vals[lb]; // Remove double-counting at start.

                if (halfWide) {
                    switch (maxType) {
                        case RECTANGULAR:
                            vals[0] = 0.5 * width * fX[0];
                            vals[numPoints - 2] += (0.5 * width * fX[numPoints - 2]);
                            break;
                        default:
                            vals[0] = 0.25 * width * fX[0];
                            vals[1] += 0.25 * width * fX[1];
                            vals[numPoints - 2] += (0.25 * width * fX[numPoints - 2]);
                            vals[numPoints - 1] = (0.25 * width * fX[numPoints - 1]);
                            break;
                    }
                }
                break;
        }

        return vals;
    }

    /**
     * Generic caller for 1D integration schemes given an IntegrationType.
     *
     * @param data To integrate
     * @param side Integrate from side
     * @param type Scheme to use
     * @return Numeric integral
     */
    public static double integrateData(DataSet data, IntegrationSide side, IntegrationType type) {
        switch (type) {
            case RECTANGULAR:
                return rectangular(data, side);
            case TRAPEZOIDAL:
                return trapezoidal(data, side);
            case SIMPSONS:
                return simpsons(data, side);
            case BOOLE:
                return booles(data, side);
            default:
                logger.warning(String.format(" Integration type %s not recognized! Defaulting to Simpson's integration", type));
                return simpsons(data, side);
        }
    }

    /**
     * Numerically integrates a data set using rectangular integration. Not
     * recommended; preferably use at least trapezoidal integration.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @return Area of integral
     */
    public static double rectangular(DataSet data, IntegrationSide side) {
        double area = 0;
        int lb = 0;
        int ub = data.numPoints() - 1;
        if (data.halfWidthEnds()) {
            area = rectangularEnds(data, side);
            ++lb;
            --ub;
        }
        area += rectangular(data, side, lb, ub);
        return area;
    }

    /**
     * Numerically integrates a data set, in bounds lb-ub inclusive, using
     * rectangular integration. Not recommended; preferably use at least
     * trapezoidal integration.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @param lb   First index to integrate over
     * @param ub   Last index to integrate over
     * @return Area of integral
     */
    public static double rectangular(DataSet data, IntegrationSide side, int lb, int ub) {
        double area = 0;
        double width = data.binWidth();
        double[] points = data.getAllFxPoints();
        assert ub > lb;
        assert ub < points.length;

        switch (side) {
            case RIGHT:
                for (int i = ub; i > lb; i--) {
                    area += (width * points[i]);
                }
                break;
            case LEFT:
            default:
                for (int i = lb; i < ub; i++) {
                    area += width * points[i];
                }
                break;
        }
        return area;
    }

    /**
     * Treats half-width bins at the ends of a DataSet using rectangular
     * integration. Not recommended, preferably use trapezoidal integration.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @return Estimate of area in half-width start/end bins.
     */
    public static double rectangularEnds(DataSet data, IntegrationSide side) {
        double width = 0.5 * data.binWidth();
        double area = 0;
        int npts = data.numPoints();
        switch (side) {
            case LEFT:
                area = data.getFxPoint(0) * width;
                area += (data.getFxPoint(npts - 2) * width);
                break;
            case RIGHT:
                area = data.getFxPoint(1) * width;
                area += (data.getFxPoint(npts - 1) * width);
                break;
        }
        return area;
    }

    /**
     * Numerically integrates a data set using parallelized rectangular
     * integration. Not recommended; preferably use at least trapezoidal
     * integration, and avoid parallelized versions unless necessary.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @return Area of integral
     */
    public static double rectangularParallel(DataSet data, IntegrationSide side) {
        double area = 0;
        int lb = 0;
        int ub = data.numPoints() - 1;
        if (data.halfWidthEnds()) {
            area = rectangularEnds(data, side);
            ++lb;
            --ub;
        }
        area += rectangularParallel(data, side, lb, ub);
        return area;
    }

    /**
     * Numerically integrates a data set, in bounds lb-ub inclusive, using
     * rectangular integration. Not recommended; preferably use at least
     * trapezoidal integration. Also, prefer parallelized versions unless
     * necessary.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @param lb   First index to integrate over
     * @param ub   Last index to integrate over
     * @return Area of integral
     */
    public static double rectangularParallel(DataSet data, IntegrationSide side, int lb, int ub) {
        double width = data.binWidth();
        double[] points = data.getAllFxPoints();
        assert ub > lb;
        assert ub < points.length;

        if (side == IntegrationSide.RIGHT) {
            ++ub;
            ++lb;
        }
        return IntStream.range(lb, ub).parallel().mapToDouble((int i) -> {
            return points[i] * width;
        }).sum();
    }

    /**
     * Numerically integrates a data set using Simpson's rule.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @return Area of integral
     */
    public static double simpsons(DataSet data, IntegrationSide side) {
        double area = 0;
        int lb = 0;
        int ub = data.numPoints() - 1;
        if (data.halfWidthEnds()) {
            area = trapezoidalEnds(data, side);
            lb++;
            ub--;
        }
        area += simpsons(data, side, lb, ub);
        return area;
    }

    /**
     * Numerically integrates a data set, in bounds lb-ub inclusive, using
     * Simpson's rule.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @param lb   First index to integrate over
     * @param ub   Last index to integrate over
     * @return Area of integral
     */
    public static double simpsons(DataSet data, IntegrationSide side, int lb, int ub) {
        double area = 0;
        double width = data.binWidth();
        double[] points = data.getAllFxPoints();
        int increment = 2;

        int nBins = (ub - lb) / increment;
        int lowerNeglected;
        int upperNeglected;

        switch (side) {
            case RIGHT:
                for (int i = ub; i > (lb + increment - 1); i -= increment) {
                    area += points[i] + (4 * points[i - 1]) + points[i - 2];
                }
                lowerNeglected = lb;
                upperNeglected = ub - (increment * nBins);
                break;
            case LEFT:
            default:
                for (int i = lb; i < (ub - increment + 1); i += increment) {
                    area += points[i] + (4 * points[i + 1]) + points[i + 2];
                }
                lowerNeglected = lb + (increment * nBins);
                upperNeglected = ub;
                break;
        }
        area *= ONE_THIRD;
        area *= width;

        area += finishIntegration(data, side, lowerNeglected, upperNeglected, IntegrationType.SIMPSONS);
        return area;
    }

    /**
     * Numerically integrates a data set using Boole's rule. The sequential
     * version is preferred unless necessary.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @return Area of integral
     */
    public static double simpsonsParallel(DataSet data, IntegrationSide side) {
        double area = 0;
        int lb = 0;
        int ub = data.numPoints() - 1;
        if (data.halfWidthEnds()) {
            area = trapezoidalEnds(data, side);
            lb++;
            ub--;
        }
        area += simpsonsParallel(data, side, lb, ub);
        return area;
    }

    /**
     * Numerically integrates a data set, in bounds lb-ub inclusive, using
     * Simpson's rule. The sequential version is preferred unless necessary.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @param lb   First index to integrate over
     * @param ub   Last index to integrate over
     * @return Area of integral
     */
    public static double simpsonsParallel(DataSet data, IntegrationSide side, int lb, int ub) {
        double area = 0;
        double width = data.binWidth();
        double[] points = data.getAllFxPoints();
        int increment = 2;

        int nBins = (ub - lb) / increment;
        int lowerNeglected;
        int upperNeglected;

        switch (side) {
            case RIGHT:
                lowerNeglected = lb;
                upperNeglected = ub - (increment * nBins);
                area = IntStream.range(0, nBins).parallel().mapToDouble((int i) -> {
                    int fromPoint = ub - (increment * i);
                    return (points[fromPoint - 2] + 4 * points[fromPoint - 1] + points[fromPoint]);
                }).sum();
                break;
            case LEFT:
            default:
                area = IntStream.range(0, nBins).parallel().mapToDouble((int i) -> {
                    int fromPoint = lb + (increment * i);
                    return (points[fromPoint] + 4 * points[fromPoint + 1] + points[fromPoint + 2]);
                }).sum();
                lowerNeglected = lb + (increment * nBins);
                upperNeglected = ub;
                break;
        }
        area *= ONE_THIRD;
        area *= width;

        area += finishIntegration(data, side, lowerNeglected, upperNeglected, IntegrationType.SIMPSONS);
        return area;
    }

    /**
     * Numerically integrates a data set using trapezoidal integration. For most
     * data sets, Simpson's rule or Boole's rule will out-perform trapezoidal
     * integration.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @return Area of integral
     */
    public static double trapezoidal(DataSet data, IntegrationSide side) {
        double area = 0;
        int lb = 0;
        int ub = data.numPoints() - 1;
        if (data.halfWidthEnds()) {
            area = trapezoidalEnds(data, side);
            lb++;
            ub--;
        }
        area += trapezoidal(data, side, lb, ub);
        return area;
    }

    /**
     * Numerically integrates a data set, in bounds lb-ub inclusive, using
     * trapezoidal integration. For most data sets, Simpson's rule or Boole's
     * rule will out-perform trapezoidal integration.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @param lb   First index to integrate over
     * @param ub   Last index to integrate over
     * @return Area of integral
     */
    public static double trapezoidal(DataSet data, IntegrationSide side, int lb, int ub) {
        double width = data.binWidth();
        double[] points = data.getAllFxPoints();

        double area = 0.5 * points[lb];
        area += 0.5 * points[ub];
        for (int i = lb + 1; i < ub; i++) {
            area += points[i];
        }
        area *= width;
        return area;
    }

    /**
     * Treats half-width bins at the ends of a DataSet using trapezoidal
     * integration.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @return Estimate of area in half-width start/end bins.
     */
    public static double trapezoidalEnds(DataSet data, IntegrationSide side) {
        double width = 0.5 * data.binWidth();
        int nPts = data.numPoints();
        double area = data.getFxPoint(0) + data.getFxPoint(1) + data.getFxPoint(nPts - 2) + data.getFxPoint(nPts - 1);
        area *= (0.5 * width);
        return area;
    }

    /**
     * Numerically integrates a data set using trapezoidal integration. For most
     * data sets, Simpson's rule or Boole's rule will out-perform trapezoidal
     * integration. Prefer use of the sequential version unless necessary.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @return Area of integral
     */
    public static double trapezoidalParallel(DataSet data, IntegrationSide side) {
        double area = 0;
        int lb = 0;
        int ub = data.numPoints() - 1;
        if (data.halfWidthEnds()) {
            area = trapezoidalEnds(data, side);
            lb++;
            ub--;
        }
        area += trapezoidalParallel(data, side, lb, ub);
        return area;
    }

    /**
     * Numerically integrates a data set, in bounds lb-ub inclusive, using
     * trapezoidal integration. For most data sets, Simpson's rule or Boole's
     * rule will out-perform trapezoidal integration. Prefer use of the
     * sequential version unless necessary.
     *
     * @param data Data set to integrate
     * @param side Side to integrate from
     * @param lb   First index to integrate over
     * @param ub   Last index to integrate over
     * @return Area of integral
     */
    public static double trapezoidalParallel(DataSet data, IntegrationSide side, int lb, int ub) {
        double width = data.binWidth();
        double[] points = data.getAllFxPoints();

        double area = 0.5 * points[lb];
        area += 0.5 * points[ub];
        area += IntStream.range(lb + 1, ub).parallel().mapToDouble((int i) -> {
            return points[i];
        }).sum();
        area *= width;
        return area;
    }

    /**
     * Integrates the remaining points after higher-order integration rules
     * cannot evenly fit over the remaining data. For example, Boole's rule
     * requires 5 points; if a data set has 7 points, it will integrate the
     * first 5, but not handle the remaining 2; this method will, apply the
     * highest-order integration rule that the remaining points (including the
     * last one used previously) permit. In that case, it would be Simpson's
     * rule.
     *
     * @param data Finish numerical integration on
     * @param side Integration side
     * @param lb   Index of lower bound to complete integration on
     * @param ub   Index of upper bound to complete integration on
     * @param type Highest-order rule permitted to finish integration
     * @return Area of un-integrated region
     */
    private static double finishIntegration(DataSet data, IntegrationSide side, int lb, int ub, IntegrationType type) {
        int totPoints = (ub - lb);

        int perBin = type.pointsNeeded();
        int increment = perBin - 1;
        increment = Math.max(1, increment); // Needed for rectangular integration

        int nBins = totPoints / increment;
        int remainder = totPoints % increment;

        IntegrateWindow intMode;
        switch (type) {
            case BOOLE:
                intMode = Integrate1DNumeric::booles;
                break;
            case SIMPSONS:
                intMode = Integrate1DNumeric::simpsons;
                break;
            case RECTANGULAR:
                intMode = Integrate1DNumeric::rectangular;
                break;
            case TRAPEZOIDAL:
            default:
                intMode = Integrate1DNumeric::trapezoidal;
                break;
        }

        double area = 0.0;

        switch (side) {
            case RIGHT:
                for (int i = ub; i > (lb - 1 + increment); i -= increment) {
                    area += intMode.toArea(data, side, (i - increment), i);
                }
                ub -= (nBins * increment);
                break;
            case LEFT:
            default:
                for (int i = lb; i < (ub + 1 - increment); i += increment) {
                    area += intMode.toArea(data, side, i, (i + increment));
                }
                lb += (nBins * increment);
                break;
        }

        assert remainder == (ub - lb);
        switch (remainder) {
            case 0:
                break;
            case 1:
                area += trapezoidal(data, side, lb, ub);
                break;
            case 2:
                area += simpsons(data, side, lb, ub);
                break;
            case 3:
                // Alternately implement Simpson's 3/8 4-point rule.
                area += simpsons(data, side, lb, ub); // Will recursively call finishIntegration, thus getting another round of trapezoidal.
                break;
            case 4:
            default:
                throw new IllegalArgumentException("This should not be currently possible.");
        }
        return area;
    }

    /**
     * Left vs right-hand integration; left-hand integration will start from the
     * first available point, run right as far as possible, and then clean up any
     * remaining points using finishIntegration, while right-hand integration
     * will start from the last available point, run left as far as possible,
     * and then clean up any remaining points using finishIntegration.
     * <p>
     * Values: LEFT, RIGHT.
     * <p>
     * Not meaningful for trapezoidal integration; there will never be any
     * unused points, and the method is symmetrical.
     */
    public enum IntegrationSide {
        LEFT, RIGHT
    }

    /**
     * Enumeration of implemented integration methods, and the number of points
     * required by them.
     */
    public enum IntegrationType {
        RECTANGULAR(1),
        TRAPEZOIDAL(2),
        SIMPSONS(3),
        BOOLE(5);

        private final int pointsNeeded;

        IntegrationType(int points) {
            pointsNeeded = points;
        }

        public final int pointsNeeded() {
            return pointsNeeded;
        }
    }

    /**
     * Functional interface used to select an integration method (primarily for
     * finishIntegration).
     */
    @FunctionalInterface
    private interface IntegrateWindow {
        /**
         * Numerically integrates a range of x given f(x).
         *
         * @param data x and f(x) data to integrate
         * @param side Side to integrate from
         * @param lb   Lower bound of x[] to integrate
         * @param ub   Upper bound of x[] to integrate
         * @return Area of integral f(x) dx, from point lb to point ub
         */
        double toArea(DataSet data, IntegrationSide side, int lb, int ub);
    }
}
