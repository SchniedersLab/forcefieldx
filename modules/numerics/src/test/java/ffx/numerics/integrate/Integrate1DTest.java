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

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import ffx.numerics.integrate.Integrate1DNumeric.IntegrationSide;
import ffx.numerics.integrate.Integrate1DNumeric.IntegrationType;

import static ffx.numerics.integrate.Integrate1DNumeric.IntegrationSide.*;
import java.util.logging.Logger;

/**
 * The IntegrationTest is a JUnit test for the Integration program that ensures
 * that the integrals are calculated correctly. This test is run using known
 * integrals calculated with the equation y=10sin(6x)-7cos(5x)+11sin(8x).
 *
 * @author ceoconnell
 */
public class Integrate1DTest {
    
    private final static int NUM_INTEGRATION_TYPES = 8; // Left/right, rect/trap/simp/boole.
    private final static Logger logger = Logger.getLogger(Integrate1DTest.class.getName());
    
    /**
     * Basic test for polynomial functions; checks to ensure that each integration
     * method is behaving as expected, and returns the exact value (to within
     * machine precision) when the method should be exact. For example, Simpson's
     * method fits a quadratic curve to 3 points, and should return the exact
     * integral for polynomials of degree 2 or less.
     */
    @Test
    public void polynomialTest() {
        logger.info(" Testing integration methods on polynomials with known integrals.");
        
        /**
         * This sets up 9 evenly spaced points from 0-1, which should be 0, 0.125,
         * 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, and 1.0. Having very few points
         * of integration should stress the integration methods to the limit...
         * not that it appears to be causing any problems for Simpson's method
         * on the third-order integral.
         * 
         * 9 points also helps because neither Simpson's nor Boole's rule should
         * be truncated and require lower-order integration at an end. If we
         * add Simpson's 3/8 rule, 13 points may be preferable.
         */
        double[] points = Integrate1DNumeric.generateXPoints(0.0, 1.0, 9, false);
        
        /**
         * These are the coefficients for polynomials of order 0 to 5, in reverse
         * order; thus the second-order polynomial is x^2 - 4x + 1.5.
         */
        double[] zeroOrder = {1.0}; // f(x) = 1
        double[] firstOrder = {2.0, 1.0}; // f(x) = x+2
        double[] secondOrder = {1.5, -4.0, 1.0}; // f(x) = x^2 - 4x + 1.5
        double[] thirdOrder = {2.0, -1.0, -3.0, 1.0}; // f(x) = x^3 - 3x^2 - x + 2
        double[] fourthOrder = {1, -4.0, -6.0, 4.0, 1.0}; // f(x) = x^4 + 4x^3 - 6x^2 - 4x + 1
        double[] fifthOrder = {2.0, 10.0, -18.0, 8.0, -5.0, 1.0}; // f(x) = x^5 - 5x^4 + 8x^3 - 18x^2 + 10x + 2
        
        /**
         * Hand-calculated integrals from 0-1, to test that the analyticalIntegral
         * function is indeed working for PolynomialCurve.
         */
        double[] trueVals = {1.0, 2.5, (-1.0 / 6.0), 0.75, -1.8, (13.0/6.0)};
        
        // Accumulate the coefficients into a 2-D array.
        double[][] coeffs = {zeroOrder, firstOrder, secondOrder, thirdOrder, fourthOrder, fifthOrder};
        
        /**
         * A bit of streaming/mapping logic here; coeffs[][] is streamed to a
         * Stream of double[]. For each of these double[], a PolynomialCurve
         * is constructed using the points array (points along x), explicitly
         * setting halfWidthBins to false, and then the streamed double[] of
         * coefficients. That transforms (maps) a Stream of double[] to a Stream
         * of PolynomialCurve, which is then collected into a list.
         * 
         * Each PolynomialCurve is:
         * A DataSet, describing a set of points f(x), starting at lb (lower bound),
         * ending at ub (upper bound), with even spacing with the optional exception
         * of half-width bins. The original set of points x[] is not guaranteed
         * to be stored, but can be reconstructed by knowing the number of points,
         * the lower bound, the upper bound, and if it used half-width start/end
         * bins.
         * 
         * A FunctionDataCurve, describing a DataSet generated from some function
         * f(x) which can be analytically integrated; this adds methods for
         * analytical integration and evaluating f(x) at any arbitrary x.
         * 
         * A PolynomialCurve, a concrete class extending FunctionDataCurve, where
         * f(x) is an arbitrary-order polynomial function.
         * 
         * This is very standard inheritance; an interface (DataSet) to describe
         * what you want something to do, an abstract class (FunctionDataCurve)
         * to provide some basic implementations, and, in this case, some 
         * additional methods for that particular subset of its interface,
         * and finally a concrete class (PolynomialCurve) which defines as little
         * as possible.
         * 
         * The whole point of this, then, is that the integration methods can
         * use any arbitrary DataSet, such as a DoublesDataSet taken straight
         * from OSRW, and can be thoroughly tested using any arbitrary 
         * FunctionDataCurve, so that we're confident they are working as intended
         * when applied to DataSets where we don't know the true answer.
         */
        
        /**
         * Create the polynomial curves, using points[] as x, disabling half-width
         * bins, and using each double[] of coefficients in turn. The values of
         * f(x) are automatically generated by the PolynomialCurve constructor.
         */
        List<FunctionDataCurve> polynomials = Arrays.stream(coeffs).map((double[] coeff) -> {
            return new PolynomialCurve(points, false, coeff);
        }).collect(Collectors.toList());
        
        for (int i = 0; i < polynomials.size(); i++) {
            FunctionDataCurve pn = polynomials.get(i);
            
            // Get the true value, and analyticalIntegral over the range.
            double trueVal = trueVals[i];
            double analytical = pn.analyticalIntegral();
            
            double rLeft = Integrate1DNumeric.rectangular(pn, LEFT);
            double rRight = Integrate1DNumeric.rectangular(pn, RIGHT);
            
            double tLeft = Integrate1DNumeric.trapezoidal(pn, LEFT);
            double tRight = Integrate1DNumeric.trapezoidal(pn, RIGHT);
            
            double sLeft = Integrate1DNumeric.simpsons(pn, LEFT);
            double sRight = Integrate1DNumeric.simpsons(pn, RIGHT);
            
            double bLeft = Integrate1DNumeric.booles(pn, LEFT);
            double bRight = Integrate1DNumeric.booles(pn, RIGHT);
            
            StringBuilder sb = new StringBuilder();
            
            sb.append(String.format(" Integrals for polynomial of degree %d\n", i));
            sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Exact", trueVal, "Analytical", analytical));
            sb.append(" Numerical integration errors:\n");
            sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Left rectangular", rLeft - trueVal, "Right rectangular", rRight - trueVal));
            sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Left trapezoidal", tLeft - trueVal, "Right trapezoidal", tRight - trueVal));
            sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Left Simpson's", sLeft - trueVal, "Right Simpson's", sRight - trueVal));
            sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Left Boole's", bLeft - trueVal, "Right Boole's", bRight - trueVal));
            sb.append("\n");
            logger.info(sb.toString());
            
            assertToUlp(trueVal, 10.0, analytical);
            if (i == 0) {
                assertToUlp(trueVal, 500.0, rLeft, rRight);
            }
            if (i <= 1) {
                assertToUlp(trueVal, 500.0, tLeft, tRight);
            }
            if (i <= 2) {
                assertToUlp(trueVal, 500.0, sLeft, sRight);
            }
            if (i <= 4) {
                assertToUlp(trueVal, 500.0, bLeft, bRight);
            }
        }
    }
    
    /**
     * A more difficult test on polynomials, with just five points and larger
     * coefficients; intended to test stability in extreme cases.
     */
    @Test
    public void polynomialGrinderTest () {
        logger.info(" Testing integration methods on polynomials with large coefficients.");
        double[] points = Integrate1DNumeric.generateXPoints(0.0, 2.0, 5, false);
        
        double[] zeroOrder = {-2.0};
        double[] firstOrder = {6.0, -14.0};
        double[] secondOrder = {4.5, -5.0, 3.0};
        double[] thirdOrder = {-3.0, 5.0, -2.0, 12.0};
        double[] fourthOrder = {12, -4.5, 8.0, -5.0, 4.0};
        double[] fifthOrder = {-4.0, 10.0, -18.0, 8.0, -5.0, 8.0};
        
        double[][] coeffs = {zeroOrder, firstOrder, secondOrder, thirdOrder, fourthOrder, fifthOrder};
        
        List<FunctionDataCurve> polynomials = Arrays.stream(coeffs).map((double[] coeff) -> {
            return new PolynomialCurve(points, false, coeff);
        }).collect(Collectors.toList());
        
        for (int i = 0; i < polynomials.size(); i++) {
            FunctionDataCurve pn = polynomials.get(i);
            
            // Get the analyticalIntegral over the range.
            double analytical = pn.analyticalIntegral();
            
            double rLeft = Integrate1DNumeric.rectangular(pn, LEFT);
            double rRight = Integrate1DNumeric.rectangular(pn, RIGHT);
            
            double tLeft = Integrate1DNumeric.trapezoidal(pn, LEFT);
            double tRight = Integrate1DNumeric.trapezoidal(pn, RIGHT);
            
            double sLeft = Integrate1DNumeric.simpsons(pn, LEFT);
            double sRight = Integrate1DNumeric.simpsons(pn, RIGHT);
            
            double bLeft = Integrate1DNumeric.booles(pn, LEFT);
            double bRight = Integrate1DNumeric.booles(pn, RIGHT);
            
            StringBuilder sb = new StringBuilder();
            
            sb.append(String.format(" Integrals for polynomial of degree %d\n", i));
            sb.append(String.format(" %-18s %9.3g\n", "Analytical", analytical));
            sb.append(" Numerical integration errors:\n");
            sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Left rectangular", rLeft - analytical, "Right rectangular", rRight - analytical));
            sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Left trapezoidal", tLeft - analytical, "Right trapezoidal", tRight - analytical));
            sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Left Simpson's", sLeft - analytical, "Right Simpson's", sRight - analytical));
            sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Left Boole's", bLeft - analytical, "Right Boole's", bRight - analytical));
            sb.append("\n");
            
            logger.info(sb.toString());
            
            if (i == 0) {
                assertToUlp(analytical, 500.0, rLeft, rRight);
            }
            if (i <= 1) {
                assertToUlp(analytical, 500.0, tLeft, tRight);
            }
            if (i <= 2) {
                assertToUlp(analytical, 500.0, sLeft, sRight);
            }
            if (i <= 4) {
                assertToUlp(analytical, 500.0, bLeft, bRight);
            }
        }
        
    }
    
    @Test
    public void testSinCosine() {
        double[] points = Integrate1DNumeric.generateXPoints(0, 1, 201, false);
        boolean verb = false;
        // Test sine with full-width ends.
        logger.info("\n Sin wave test without half-width bins");
        sinTest(points,  false, false, verb);
        logger.info("\n Cosine wave test without half-width bins");
        sinTest(points, false, true, verb);
        
        points = Integrate1DNumeric.generateXPoints(0, 1, 202, true);
        logger.info("\n Sin wave test with half-width bins");
        sinTest(points, true, false, verb);
        logger.info("\n Cosine wave test with half-width bins");
        sinTest(points, true, true, verb);
        logger.info("\n");
    }
    
    /**
     * Common code for testing sin/cosine waves; effectively a sub-method for 
     * the testSinCosine test.
     * @param points
     * @param halvedEnds
     * @param cosine If cosine, else sine wave
     * @param verbose 
     */
    public void sinTest(double[] points, boolean halvedEnds, boolean cosine, boolean verbose) {
        int failRleft = 0;
        int failRright = 0;
        int failTleft = 0;
        int failTright = 0;
        int failSleft = 0;
        int failSright = 0;
        int failBleft = 0;
        int failBright = 0;
        double maxDelta = 0.05;
        
        for (int j = 1; j <= 500; j++) {
            FunctionDataCurve wave = cosine ? new CosineWave(points, halvedEnds, j, j) : new SinWave(points, halvedEnds, j, j);
            
            // Get the analyticalIntegral over the range.
            double analytical = wave.analyticalIntegral();
            
            double rLeft = Integrate1DNumeric.rectangular(wave, LEFT);
            double erLeft = analytical - rLeft;
            double rRight = Integrate1DNumeric.rectangular(wave, RIGHT);
            double erRight = analytical - rRight;
            
            double tLeft = Integrate1DNumeric.trapezoidal(wave, LEFT);
            double etLeft = analytical - tLeft;
            double tRight = Integrate1DNumeric.trapezoidal(wave, RIGHT);
            double etRight = analytical - tRight;
            
            double sLeft = Integrate1DNumeric.simpsons(wave, LEFT);
            double esLeft = analytical - sLeft;
            double sRight = Integrate1DNumeric.simpsons(wave, RIGHT);
            double esRight = analytical - sRight;
            
            double bLeft = Integrate1DNumeric.booles(wave, LEFT);
            double ebLeft = analytical - bLeft;
            double bRight = Integrate1DNumeric.booles(wave, RIGHT);
            double ebRight = analytical - bRight;
            if (verbose) {
                StringBuilder sb = new StringBuilder();
                
                sb.append(String.format(" Integrals for sine wave j*sin(j) of j %d\n", j));
                sb.append(String.format(" %-18s %9.3g\n", "Analytical", analytical));
                sb.append(" Numerical integration errors:\n");
                sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Left rectangular", erLeft, "Right rectangular", erRight));
                sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Left trapezoidal", etLeft, "Right trapezoidal", etRight));
                sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Left Simpson's", esLeft, "Right Simpson's", esRight));
                sb.append(String.format(" %-18s %9.3g  %-18s %9.3g\n", "Left Boole's", ebLeft, "Right Boole's", ebRight));
                sb.append("\n");
                
                logger.info(sb.toString());
            }
            
            if (failRleft == 0 && Math.abs(erLeft) > maxDelta) {
                failRleft = j;
            }
            if (failRright == 0 && Math.abs(erRight) > maxDelta) {
                failRright = j;
            }
            
            if (failTleft == 0 && Math.abs(etLeft) > maxDelta) {
                failTleft = j;
            }
            if (failTright == 0 && Math.abs(etRight) > maxDelta) {
                failTright = j;
            }
            
            if (failSleft == 0 && Math.abs(esLeft) > maxDelta) {
                failSleft = j;
            }
            if (failSright == 0 && Math.abs(esRight) > maxDelta) {
                failSright = j;
            }
            
            if (failBleft == 0 && Math.abs(ebLeft) > maxDelta) {
                failBleft = j;
            }
            if (failBright == 0 && Math.abs(ebRight) > maxDelta) {
                failBright = j;
            }
            
        }
        
        StringBuilder sb1 = new StringBuilder(String.format(" Error exceeded delta %8.3f at iterations:\n", maxDelta));
        sb1.append(String.format(" %-18s %9d  %-18s %9d\n", "Left rectangular", failRleft, "Right rectangular", failRright));
        sb1.append(String.format(" %-18s %9d  %-18s %9d\n", "Left trapezoidal", failTleft, "Right trapezoidal", failTright));
        sb1.append(String.format(" %-18s %9d  %-18s %9d\n", "Left Simpson's", failSleft, "Right Simpson's", failSright));
        sb1.append(String.format(" %-18s %9d  %-18s %9d\n", "Left Boole's", failBleft, "Right Boole's", failBright));
        logger.info(sb1.toString());
        
        StringBuilder sb = new StringBuilder(" integration failed for ");
        if (halvedEnds) {
            if (cosine) {
                sb.append("cosine wave, with half-width end bins, of frequency ");
            } else {
                sb.append("sine wave, with half-width end bins, of frequency ");
            }
            String midMessage = sb.toString();
            assertTrue(String.format("Rectangular%s%d", midMessage, 12), failRleft > 12 && failRright > 12);
            assertTrue(String.format("Trapezoidal%s%d", midMessage, 100), failTleft > 100 && failTright > 100);
            assertTrue(String.format("Simpson's%s%d", midMessage, 150), failSleft > 150 && failSright > 150);
            assertTrue(String.format("Boole's%s%d", midMessage, 150), failBleft > 150 && failBright > 150);
            
        } else {
            if (cosine) {
                sb.append("cosine wave of frequency ");
            } else {
                sb.append("sine wave of frequency ");
            }
            String midMessage = sb.toString();
            assertTrue(String.format("Rectangular%s%d", midMessage, 12), failRleft > 12 && failRright > 12);
            assertTrue(String.format("Trapezoidal%s%d", midMessage, 100), failTleft > 100 && failTright > 100);
            assertTrue(String.format("Simpson's%s%d", midMessage, 250), failSleft > 250 && failSright > 250);
            assertTrue(String.format("Boole's%s%d", midMessage, 240), failBleft > 240 && failBright > 240);
        }
    }
    
    /**
     * Ensures that integration performs the same on a DoublesDataSet constructed
     * from a function as directly from the function; uses both DoublesDataSet
     * constructors.
     */
    @Test
    public void castToDoubleSetTest() {
        logger.info(" Testing casting of a function to a DoublesDataSet");
        
        double[] points = Integrate1DNumeric.generateXPoints(0.0, 1.0, 11, false);
        
        double[] zeroOrder = {1.0}; // f(x) = 1
        double[] firstOrder = {2.0, 1.0}; // f(x) = x+2
        double[] secondOrder = {1.5, -4.0, 1.0}; // f(x) = x^2 - 4x + 1.5
        double[] thirdOrder = {2.0, -1.0, -3.0, 1.0}; // f(x) = x^3 - 3x^2 - x + 2
        double[] fourthOrder = {1, -4.0, -6.0, 4.0, 1.0}; // f(x) = x^4 + 4x^3 - 6x^2 - 4x + 1
        double[] fifthOrder = {2.0, 10.0, -18.0, 8.0, -5.0, 1.0}; // f(x) = x^5 - 5x^4 + 8x^3 - 18x^2 + 10x + 2
        
        // Accumulate the coefficients into a 2-D array.
        double[][] coeffs = {zeroOrder, firstOrder, secondOrder, thirdOrder, fourthOrder, fifthOrder};
        List<FunctionDataCurve> polynomials = Arrays.stream(coeffs).map((double[] coeff) -> {
            return new PolynomialCurve(points, false, coeff);
        }).collect(Collectors.toList());
        
        for (FunctionDataCurve pn : polynomials) {
            DataSet dpn = new DoublesDataSet(pn);
            double[] ypts = pn.getAllFxPoints();
            DataSet abInitio = new DoublesDataSet(points, ypts, false);
            
            for (int j = 0; j < NUM_INTEGRATION_TYPES; j++) {
                IntegrationResult rpn = new IntegrationResult(pn, j);
                IntegrationResult rdpn = new IntegrationResult(dpn, j);
                IntegrationResult raipn = new IntegrationResult(abInitio, j);
                
                double pnVal = rpn.getValue();
                double dpnVal = rdpn.getValue();
                double aipnVal = raipn.getValue();
                
                String message = String.format(" Casted function:\n%s\nfunction "
                        + "result %9.3g, casted result %9.3g, from points %9.3g, "
                        + "difference to casted %9.3g, difference to from-points %9.3g"
                        + "\nintegration method %s", pn.toString(), pnVal, dpnVal, 
                        aipnVal, (pnVal - dpnVal), (pnVal - aipnVal), rpn.toString());
                assertToUlp(message, pnVal, 10.0, dpnVal, aipnVal);
            }
        }
        logger.info(" Polynomial casting tests without halved sides complete.");
        
        double[] points2 = Integrate1DNumeric.generateXPoints(0.0, 1.0, 12, true);
        polynomials = Arrays.stream(coeffs).map((double[] coeff) -> {
            return new PolynomialCurve(points2, true, coeff);
        }).collect(Collectors.toList());
        
        for (FunctionDataCurve pn : polynomials) {
            DataSet dpn = new DoublesDataSet(pn);
            double[] ypts = pn.getAllFxPoints();
            DataSet abInitio = new DoublesDataSet(points2, ypts, true);
            
            for (int j = 0; j < NUM_INTEGRATION_TYPES; j++) {
                IntegrationResult rpn = new IntegrationResult(pn, j);
                IntegrationResult rdpn = new IntegrationResult(dpn, j);
                IntegrationResult raipn = new IntegrationResult(abInitio, j);
                
                double pnVal = rpn.getValue();
                double dpnVal = rdpn.getValue();
                double aipnVal = raipn.getValue();
                
                String message = String.format(" Casted function:\n%s\nfunction "
                        + "result %9.3g, casted result %9.3g, from points %9.3g, "
                        + "difference to casted %9.3g, difference to from-points %9.3g"
                        + "\nintegration method %s", pn.toString(), pnVal, dpnVal, 
                        aipnVal, (pnVal - dpnVal), (pnVal - aipnVal), rpn.toString());
                assertToUlp(message, pnVal, 10.0, dpnVal, aipnVal);
            }
        }
        logger.info(" Polynomial casting tests with halved sides complete.");
        
        double[] sinePoints = Integrate1DNumeric.generateXPoints(0, 4.0, 101, false);
        FunctionDataCurve sinWave = new SinWave(sinePoints, false, 20, 30);
        DataSet dsw = new DoublesDataSet(sinWave);
        double[] ypts = sinWave.getAllFxPoints();
        DataSet abInitio = new DoublesDataSet(sinePoints, ypts, false);
        
        for (int j = 0; j < NUM_INTEGRATION_TYPES; j++) {
            IntegrationResult rsw = new IntegrationResult(sinWave, j);
            IntegrationResult rdsw = new IntegrationResult(dsw, j);
            IntegrationResult raisw = new IntegrationResult(abInitio, j);
            
            double swVal = rsw.getValue();
            double dswVal = rdsw.getValue();
            double aiswVal = raisw.getValue();
                
            String message = String.format(" Casted function:\n%s\nfunction "
                    + "result %9.3g, casted result %9.3g, from points %9.3g, "
                    + "difference to casted %9.3g, difference to from-points %9.3g"
                    + "\nintegration method %s", sinWave.toString(), swVal, dswVal, 
                    aiswVal, (swVal - dswVal), (swVal - aiswVal), rsw.toString());
            assertToUlp(message, swVal, 10.0, dswVal, aiswVal);
        }
        logger.info(" Sine-wave casting tests with halved sides complete.");
        
        double[] sinePoints2 = Integrate1DNumeric.generateXPoints(0, 4.0, 102, true);
        sinWave = new SinWave(sinePoints2, true, 20, 30);
        dsw = new DoublesDataSet(sinWave);
        ypts = sinWave.getAllFxPoints();
        abInitio = new DoublesDataSet(sinePoints2, ypts, true);
        
        for (int j = 0; j < NUM_INTEGRATION_TYPES; j++) {
            IntegrationResult rsw = new IntegrationResult(sinWave, j);
            IntegrationResult rdsw = new IntegrationResult(dsw, j);
            IntegrationResult raisw = new IntegrationResult(abInitio, j);
            
            double swVal = rsw.getValue();
            double dswVal = rdsw.getValue();
            double aiswVal = raisw.getValue();
                
            String message = String.format(" Casted function:\n%s\nfunction "
                    + "result %9.3g, casted result %9.3g, from points %9.3g, "
                    + "difference to casted %9.3g, difference to from-points %9.3g"
                    + "\nintegration method %s", sinWave.toString(), swVal, dswVal, 
                    aiswVal, (swVal - dswVal), (swVal - aiswVal), rsw.toString());
            assertToUlp(message, swVal, 10.0, dswVal, aiswVal);
        }
        logger.info(" Sine-wave casting tests with halved sides complete.");
        logger.info(" Casting tests complete.");
    }
    
    /**
     * Assert that doubles are equal to within a multiplier of ulp (machine precision).
     * @param trueVal True answer
     * @param ulpMult Multiple of ulp to use
     * @param values Values to check
     */
    private static void assertToUlp(double trueVal, double ulpMult, double... values) {
        double ulp = Math.ulp(trueVal) * ulpMult;
        for (double val : values) {
            assertEquals(trueVal, val, ulp);
        }
    }
    
    /**
     * Assert that doubles are equal to within a multiplier of ulp (machine precision).
     * @param message To print if assertion fails
     * @param trueVal True answer
     * @param ulpMult Multiple of ulp to use
     * @param values Values to check
     */
    private static void assertToUlp(String message, double trueVal, double ulpMult, double... values) {
        double ulp = Math.ulp(trueVal) * ulpMult;
        for (double val : values) {
            assertEquals(message, trueVal, val, ulp);
        }
    }
    
    /**
     * Tests the parallelized versions of the integration methods. Overall,
     * these are not recommended for production use; it would require an enormous
     * number of points for numerical integration of a 1-D function to be at
     * all significant for timing, and the parallelized versions are probably
     * not quite as CPU-efficient.
     */
    @Test
    public void parallelTest() {
        double[] pts = Integrate1DNumeric.generateXPoints(0, 2.0, 92, false);
        FunctionDataCurve tcurve = new SinWave(pts, false, 60, 60);
        for (int i = 0; i < NUM_INTEGRATION_TYPES; i++) {
            IntegrationResult seqResult = new IntegrationResult(tcurve, i, false);
            IntegrationResult parResult = new IntegrationResult(tcurve, i, true);
            double seqVal = seqResult.getValue();
            double parVal = parResult.getValue();
            String message = String.format("Parallel value %9.3g did not match sequential value %9.3g, error %9.3g", parVal, seqVal, (parVal - seqVal));
            assertToUlp(message, seqVal, 80.0, parVal);
        }
    }
    
    /**
     * Intended to be a shorthand way of performing all the standard 
     * integrations, instead of tediously coding in separate tests for all
     * four integration types in both directions; auto-assigns several fields
     * based on the provided index.
     */
    private class IntegrationResult {
        private final IntegrationType type;
        private final IntegrationSide side;
        private final boolean halvedSides;
        private final DataSet data;
        private final double integratedValue;
        private final boolean parallel;
        
        /**
         * Basic constructor, automatically setting direction and integration
         * method based on index with a useful toString.
         * 
         * Even index values produce a left-hand integral, odd a right-hand 
         * integral.
         * 
         * Integration type is as follows, 0-7, after modulo(8)
         * 0-1: Rectangular
         * 2-3: Trapezoidal
         * 4-5: Simpson's rule
         * 6-7: Boole's rule
         * 
         * @param data Data to numerically integrate
         * @param index Automatically set method and direction
         */
        public IntegrationResult(DataSet data, int index) {
            this(data, index, false);
        }
        
        public IntegrationResult(DataSet data, int index, boolean parallel) {
            switch ((index / 2) % 4) {
                case 0:
                    type = IntegrationType.RECTANGULAR;
                    break;
                case 1:
                    type = IntegrationType.TRAPEZOIDAL;
                    break;
                case 2:
                    type = IntegrationType.SIMPSONS;
                    break;
                case 3:
                    type = IntegrationType.BOOLE;
                    break;
                default:
                    throw new IllegalArgumentException(" How did ((index / 2) % 4) ever come out to not be 0-3?");
            }
            side = (index % 2 == 0) ? LEFT : RIGHT;
            halvedSides = data.halfWidthEnds();
            this.data = data;
            this.parallel = parallel;
            
            if (parallel) {
                switch(type) {
                    case RECTANGULAR:
                        integratedValue = Integrate1DNumeric.rectangularParallel(data, side);
                        break;
                    case TRAPEZOIDAL:
                        integratedValue = Integrate1DNumeric.trapezoidalParallel(data, side);
                        break;
                    case SIMPSONS:
                        integratedValue = Integrate1DNumeric.simpsonsParallel(data, side);
                        break;
                    case BOOLE:
                        integratedValue = Integrate1DNumeric.boolesParallel(data, side);
                        break;
                    default:
                        throw new IllegalArgumentException(" How did this end up with an integration type not rectangular, trapezoidal, Simpson's, or Boole's?");
                }
                
            } else {
                switch(type) {
                    case RECTANGULAR:
                        integratedValue = Integrate1DNumeric.rectangular(data, side);
                        break;
                    case TRAPEZOIDAL:
                        integratedValue = Integrate1DNumeric.trapezoidal(data, side);
                        break;
                    case SIMPSONS:
                        integratedValue = Integrate1DNumeric.simpsons(data, side);
                        break;
                    case BOOLE:
                        integratedValue = Integrate1DNumeric.booles(data, side);
                        break;
                    default:
                        throw new IllegalArgumentException(" How did this end up with an integration type not rectangular, trapezoidal, Simpson's, or Boole's?");
                }
            }
        }
        
        public IntegrationType getType() {
            return type;
        }
        
        public IntegrationSide getSide() {
            return side;
        }
        
        public DataSet getDataSet() {
            return data;
        }
        
        public boolean getHalvedSides() {
            return halvedSides;
        }
        
        public double getValue() {
            return integratedValue;
        }
        
        public boolean isParallel() {
            return parallel;
        }
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder(type.toString().toLowerCase());
            sb.append(", ").append(side.toString().toLowerCase());
            sb.append("-side integral");
            if (halvedSides) {
                sb.append(" with half-width ending bins");
            }
            sb.append(".");
            return sb.toString();
        }
    }
}
