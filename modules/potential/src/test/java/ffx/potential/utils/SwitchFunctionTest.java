/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.potential.utils;

import java.util.logging.Logger;
import static java.lang.Math.ulp;

import org.apache.commons.math3.util.FastMath;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import ffx.numerics.PowerSwitch;
import ffx.numerics.SquaredTrigSwitch;
import ffx.numerics.UnivariateSwitchingFunction;
import ffx.potential.nonbonded.MultiplicativeSwitch;

/**
 * Test the various switching functions.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class SwitchFunctionTest {
    private static final Logger logger = Logger.getLogger(SwitchFunctionTest.class.getName());
    
    // A set of standard multiples of ulp(0) and ulp(1)
    private final static double ULP_ONE_2 = 2.0 * Math.ulp(1.0);
    private final static double ULP_ONE_10 = 5.0 * ULP_ONE_2;
    private final static double ULP_ONE_100 = 50.0 * ULP_ONE_2;
    
    private final static double ULP_ZERO_2 = 2.0 * Math.ulp(0.0);
    private final static double ULP_ZERO_10 = 5.0 * ULP_ZERO_2;
    private final static double ULP_ZERO_100 = 50.0 * ULP_ZERO_2;
    
    private final static double LOOSE_TOLERANCE = 0.000001;
    private final static double MID_TOLERANCE = 1.0E-10;
    
    /**
     * Tests interpolation via the PowerSwitch class.
     */
    @Test
    public void testPowerInterpolation() {
        logger.info(" Testing default power switch");
        PowerSwitch funct = new PowerSwitch();
        standardTest(funct);
        assertEquals("Default power-switch zero bound != 0.0", 0.0, funct.getZeroBound(), ULP_ZERO_2);
        assertEquals("Default power-switch one bound != 1.0", 1.0, funct.getOneBound(), ULP_ONE_2);
        assertFalse("Power switches are not constant outside the bounds.", funct.constantOutsideBounds());
        assertFalse("Power switches are not valid outside the bounds.", funct.validOutsideBounds());
        assertEquals("Default power-switch max-zero-derivative should return 0", 0, funct.getHighestOrderZeroDerivative());
        assertTrue("Default power-switch should be equal unity with symmetric inputs", funct.symmetricToUnity());
        
        for (double x = 0; x <= 1.0; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double valAt = funct.valueAt(x);
            assertEquals(String.format("Value of default power-switch at %8.4g should be itself, was %8.4g", x, valAt), x, valAt, delta);
            double derivAt = funct.firstDerivative(x);
            assertEquals(String.format("First derivative of default power switch at %8.4g should always be 1.0, was %8.4g", x, derivAt), 1.0, derivAt, ULP_ONE_10);
            double d2 = funct.secondDerivative(x);
            assertEquals(String.format("Second derivative of default power switch at %8.4g should always be 0.0, was %8.4g", x, d2), 0.0, d2, ULP_ZERO_10);
        }
        
        logger.info(" Testing manually-constructed default power switch");
        funct = new PowerSwitch(1.0, 1.0);
        standardTest(funct);
        assertEquals("Default power-switch zero bound != 0.0", 0.0, funct.getZeroBound(), ULP_ZERO_2);
        assertEquals("Default power-switch one bound != 1.0", 1.0, funct.getOneBound(), ULP_ONE_2);
        assertFalse("Power switches are not constant outside the bounds.", funct.constantOutsideBounds());
        assertFalse("Power switches are not valid outside the bounds.", funct.validOutsideBounds());
        assertEquals("Default power-switch max-zero-derivative should return 0", 0, funct.getHighestOrderZeroDerivative());
        assertTrue("Default power-switch should be equal unity with symmetric inputs", funct.symmetricToUnity());
        
        for (double x = 0; x <= 1.0; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double valAt = funct.valueAt(x);
            assertEquals(String.format("Value of default power-switch at %8.4g should be itself, was %8.4g", x, valAt), x, valAt, delta);
            double derivAt = funct.firstDerivative(x);
            assertEquals(String.format("First derivative of default power switch at %8.4g should always be 1.0, was %8.4g", x, derivAt), 1.0, derivAt, ULP_ONE_10);
            double d2 = funct.secondDerivative(x);
            assertEquals(String.format("Second derivative of default power switch at %8.4g should always be 0.0, was %8.4g", x, d2), 0.0, d2, ULP_ZERO_10);
        }
        
        logger.info(" Testing linear power switch with doubled bounds");
        funct = new PowerSwitch(0.5, 1.0);
        standardTest(funct);
        assertEquals(String.format("Power-switch %s zero bound != 0.0", funct.toString()), 0.0, funct.getZeroBound(), ULP_ZERO_2);
        assertEquals(String.format("Power-switch %s one bound != 2.0", funct.toString()), 2.0, funct.getOneBound(), 2.0*ulp(2.0));
        assertFalse("Power switches are not constant outside the bounds.", funct.constantOutsideBounds());
        assertFalse("Power switches are not valid outside the bounds.", funct.validOutsideBounds());
        assertEquals(String.format("Power-switch %s max-zero-derivative should return 0", funct.toString()), 0, funct.getHighestOrderZeroDerivative());
        assertTrue(String.format("Power-switch %s should be equal unity with symmetric inputs", funct.toString()), funct.symmetricToUnity());
        
        for (double x = 0; x <= 1.0; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double valAt = funct.valueAt(x);
            assertEquals(String.format("Value of power-switch %s at %8.4g should be 0.5 * itself, was %8.4g", funct.toString(), x, valAt), 0.5 * x, valAt, delta);
            double derivAt = funct.firstDerivative(x);
            assertEquals(String.format("First derivative of power-switch %s at %8.4g should always be 0.5, was %8.4g", funct.toString(), x, derivAt), 0.5, derivAt, ULP_ONE_10);
            double d2 = funct.secondDerivative(x);
            assertEquals(String.format("Second derivative of power-switch %s at %8.4g should always be 0.0, was %8.4g", funct.toString(), x, d2), 0.0, d2, ULP_ZERO_10);
        }
        
        logger.info(" Testing power-2 switching function");
        funct = new PowerSwitch(1.0, 2.0);
        standardTest(funct);
        assertEquals(String.format("Power-switch %s zero bound != 0.0", funct.toString()), 0.0, funct.getZeroBound(), ULP_ZERO_2);
        assertEquals(String.format("Power-switch %s one bound != 1.0", funct.toString()), 1.0, funct.getOneBound(), ULP_ONE_2);
        assertFalse("Power switches are not constant outside the bounds.", funct.constantOutsideBounds());
        assertFalse("Power switches are not valid outside the bounds.", funct.validOutsideBounds());
        assertEquals(String.format("Power-switch %s max-zero-derivative should return 0", funct.toString()), 0, funct.getHighestOrderZeroDerivative());
        
        double beta = funct.getExponent();
        for (double x = 0; x <= 1.0; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double valAt = funct.valueAt(x);
            double trueVal = x*x;
            assertEquals(String.format("Value of power-switch %s at %8.4g should be %8.4g, was %8.4g", funct.toString(), x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = funct.firstDerivative(x);
            trueVal = beta * FastMath.pow(x, (beta-1.0));
            assertEquals(String.format("First derivative of power-switch %s at %8.4g should be %8.4g, was %8.4g", funct.toString(), x, trueVal, derivAt), trueVal, derivAt, 10.0*ulp(trueVal));
            
            //trueVal = beta * (beta - 1.0) * FastMath.pow(x, (beta-2.0));
            trueVal = 2;
            double d2 = funct.secondDerivative(x);
            assertEquals(String.format("Second derivative of power-switch %s at %8.4g should always be %8.4g, was %8.4g", funct.toString(), x, trueVal, d2), trueVal, d2, ULP_ZERO_10);
        }
        
        logger.info(" Testing power-2 switching function with double-wide bounds");
        funct = new PowerSwitch(0.5, 2.0);
        double ub = funct.getOneBound();
        standardTest(funct);
        assertEquals(String.format("Power-switch %s zero bound != 0.0", funct.toString()), 0.0, funct.getZeroBound(), ULP_ZERO_2);
        assertEquals(String.format("Power-switch %s one bound != 2.0", funct.toString()), 2.0, funct.getOneBound(), 2.0*ulp(2.0));
        assertFalse("Power switches are not constant outside the bounds.", funct.constantOutsideBounds());
        assertFalse("Power switches are not valid outside the bounds.", funct.validOutsideBounds());
        assertEquals(String.format("Power-switch %s max-zero-derivative should return 0", funct.toString()), 0, funct.getHighestOrderZeroDerivative());
        
        beta = funct.getExponent();
        for (double x = 0; x <= ub; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double valAt = funct.valueAt(x);
            double trueVal = x*x*0.25;
            assertEquals(String.format("Value of power-switch %s at %8.4g should be %8.4g, was %8.4g", funct.toString(), x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = funct.firstDerivative(x);
            trueVal = beta * 0.25 * FastMath.pow(x, (beta-1.0));
            assertEquals(String.format("First derivative of power-switch %s at %8.4g should be %8.4g, was %8.4g", funct.toString(), x, trueVal, derivAt), trueVal, derivAt, 10.0*ulp(trueVal));
            
            //trueVal = beta * (beta - 1.0) * FastMath.pow(x, (beta-2.0));
            trueVal = 0.5;
            double d2 = funct.secondDerivative(x);
            assertEquals(String.format("Second derivative of power-switch %s at %8.4g should always be %8.4g, was %8.4g", funct.toString(), x, trueVal, d2), trueVal, d2, ULP_ZERO_10);
        }
        
        logger.info(" Testing power-4 switching function");
        funct = new PowerSwitch(1.0, 4.0);
        ub = funct.getOneBound();
        standardTest(funct);
        assertEquals(String.format("Power-switch %s zero bound != 0.0", funct.toString()), 0.0, funct.getZeroBound(), ULP_ZERO_2);
        assertEquals(String.format("Power-switch %s one bound != 1.0", funct.toString()), 1.0, funct.getOneBound(), ULP_ONE_2);
        assertFalse("Power switches are not constant outside the bounds.", funct.constantOutsideBounds());
        assertFalse("Power switches are not valid outside the bounds.", funct.validOutsideBounds());
        assertEquals(String.format("Power-switch %s max-zero-derivative should return 0", funct.toString()), 0, funct.getHighestOrderZeroDerivative());
        
        beta = funct.getExponent();
        for (double x = 0; x <= ub; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double valAt = funct.valueAt(x);
            double trueVal = x*x*x*x;
            assertEquals(String.format("Value of power-switch %s at %8.4g should be %8.4g, was %8.4g", funct.toString(), x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = funct.firstDerivative(x);
            trueVal = beta * FastMath.pow(x, (beta-1.0));
            assertEquals(String.format("First derivative of power-switch %s at %8.4g should be %8.4g, was %8.4g", funct.toString(), x, trueVal, derivAt), trueVal, derivAt, 10.0*ulp(trueVal));
            
            trueVal = beta * (beta - 1.0) * FastMath.pow(x, (beta-2.0));
            double d2 = funct.secondDerivative(x);
            assertEquals(String.format("Second derivative of power-switch %s at %8.4g should always be %8.4g, was %8.4g", funct.toString(), x, trueVal, d2), trueVal, d2, ULP_ZERO_10);
        }
        
        logger.info(" Testing power-4 switching function with double-wide bounds");
        funct = new PowerSwitch(0.5, 4.0);
        ub = funct.getOneBound();
        standardTest(funct);
        assertEquals(String.format("Power-switch %s zero bound != 0.0", funct.toString()), 0.0, funct.getZeroBound(), ULP_ZERO_2);
        assertEquals(String.format("Power-switch %s one bound != 2.0", funct.toString()), 2.0, funct.getOneBound(), 2.0*ulp(2.0));
        assertFalse("Power switches are not constant outside the bounds.", funct.constantOutsideBounds());
        assertFalse("Power switches are not valid outside the bounds.", funct.validOutsideBounds());
        assertEquals(String.format("Power-switch %s max-zero-derivative should return 0", funct.toString()), 0, funct.getHighestOrderZeroDerivative());
        
        beta = funct.getExponent();
        for (double x = 0; x <= ub; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double valAt = funct.valueAt(x);
            double trueVal = 0.0625*x*x*x*x;
            assertEquals(String.format("Value of power-switch %s at %8.4g should be %8.4g, was %8.4g", funct.toString(), x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = funct.firstDerivative(x);
            trueVal = beta * 0.0625 * FastMath.pow(x, (beta-1.0));
            assertEquals(String.format("First derivative of power-switch %s at %8.4g should be %8.4g, was %8.4g", funct.toString(), x, trueVal, derivAt), trueVal, derivAt, 10.0*ulp(trueVal));
            
            trueVal = 0.0625 * beta * (beta - 1.0) * FastMath.pow(x, (beta-2.0));
            double d2 = funct.secondDerivative(x);
            assertEquals(String.format("Second derivative of power-switch %s at %8.4g should always be %8.4g, was %8.4g", funct.toString(), x, trueVal, d2), trueVal, d2, ULP_ZERO_10);
        }
        
        logger.info(" Testing square-root switching function");
        funct = new PowerSwitch(1.0, 0.5);
        ub = funct.getOneBound();
        standardTest(funct);
        assertEquals(String.format("Power-switch %s zero bound != 0.0", funct.toString()), 0.0, funct.getZeroBound(), ULP_ZERO_2);
        assertEquals(String.format("Power-switch %s one bound != 1.0", funct.toString()), 1.0, funct.getOneBound(), ULP_ONE_2);
        assertFalse("Power switches are not constant outside the bounds.", funct.constantOutsideBounds());
        assertFalse("Power switches are not valid outside the bounds.", funct.validOutsideBounds());
        assertEquals(String.format("Power-switch %s max-zero-derivative should return 0", funct.toString()), 0, funct.getHighestOrderZeroDerivative());
        
        beta = funct.getExponent();
        for (double x = 0; x <= ub; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double valAt = funct.valueAt(x);
            double trueVal = FastMath.sqrt(x);
            assertEquals(String.format("Value of power-switch %s at %8.4g should be %8.4g, was %8.4g", funct.toString(), x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = funct.firstDerivative(x);
            trueVal = beta * FastMath.pow(x, (beta-1.0));
            assertEquals(String.format("First derivative of power-switch %s at %8.4g should be %8.4g, was %8.4g", funct.toString(), x, trueVal, derivAt), trueVal, derivAt, 10.0*ulp(trueVal));
            
            trueVal = beta * (beta - 1.0) * FastMath.pow(x, (beta-2.0));
            double d2 = funct.secondDerivative(x);
            assertEquals(String.format("Second derivative of power-switch %s at %8.4g should always be %8.4g, was %8.4g", funct.toString(), x, trueVal, d2), trueVal, d2, ULP_ZERO_10);
        }
        
        logger.info(" Testing square-root switching function with double-wide bounds");
        funct = new PowerSwitch(0.5, 0.5);
        ub = funct.getOneBound();
        standardTest(funct);
        assertEquals(String.format("Power-switch %s zero bound != 0.0", funct.toString()), 0.0, funct.getZeroBound(), ULP_ZERO_2);
        assertEquals(String.format("Power-switch %s one bound != 2.0", funct.toString()), 2.0, funct.getOneBound(), 2.0*ulp(2.0));
        assertFalse("Power switches are not constant outside the bounds.", funct.constantOutsideBounds());
        assertFalse("Power switches are not valid outside the bounds.", funct.validOutsideBounds());
        assertEquals(String.format("Power-switch %s max-zero-derivative should return 0", funct.toString()), 0, funct.getHighestOrderZeroDerivative());
        
        beta = funct.getExponent();
        double sqrt05 = FastMath.sqrt(0.5);
        for (double x = 0; x <= ub; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double valAt = funct.valueAt(x);
            double trueVal = sqrt05 * FastMath.sqrt(x);
            assertEquals(String.format("Value of power-switch %s at %8.4g should be %8.4g, was %8.4g", funct.toString(), x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = funct.firstDerivative(x);
            trueVal = 0.5 * beta * FastMath.pow(0.5*x, (beta - 1.0));
            assertEquals(String.format("First derivative of power-switch %s at %8.4g should be %8.4g, was %8.4g", funct.toString(), x, trueVal, derivAt), trueVal, derivAt, 10.0*ulp(trueVal));
            
            trueVal = 0.25 * beta * (beta - 1.0) * FastMath.pow(0.5*x, (beta - 2.0));
            double d2 = funct.secondDerivative(x);
            assertEquals(String.format("Second derivative of power-switch %s at %8.4g should always be %8.4g, was %8.4g", funct.toString(), x, trueVal, d2), trueVal, d2, ULP_ZERO_10);
        }
    }
    
    /**
     * Tests multiplicative switch functions.
     */
    @Test
    public void multSwitchTest() {
        logger.info(" Testing multiplicative switch functionality");
        MultiplicativeSwitch sf = new MultiplicativeSwitch();
        standardTest(sf);
        
        assertEquals("Default multiplicative switch zero bound != 0.0", 0.0, sf.getZeroBound(), ULP_ZERO_2);
        assertEquals("Default power-switch one bound != 1.0", 1.0, sf.getOneBound(), ULP_ONE_2);
        assertFalse("Power switches are not constant outside the bounds.", sf.constantOutsideBounds());
        assertFalse("Power switches are not valid outside the bounds.", sf.validOutsideBounds());
        assertEquals("Default power-switch max-zero-derivative should return 2", 2, sf.getHighestOrderZeroDerivative());
        assertTrue("Default power-switch should be equal unity with symmetric inputs", sf.symmetricToUnity());
        
        sf = new MultiplicativeSwitch(1.0, 0.0);
        standardTest(sf);
        
        sf = new MultiplicativeSwitch(9.0, 7.2);
        standardTest(sf, LOOSE_TOLERANCE);
    }
    
    @Test
    public void trigTest() {
        logger.info(" Testing trigonometric switch functionality");
        double piOverTwo = Math.PI * 0.5;
        
        SquaredTrigSwitch sf = new SquaredTrigSwitch(false);
        standardTest(sf, MID_TOLERANCE);
        double a = piOverTwo;
        assertEquals("Default sine switch zero bound != 0.0", 0.0, sf.getZeroBound(), ULP_ZERO_2);
        assertEquals("Default sine switch one bound != 1.0", 1.0, sf.getOneBound(), ULP_ONE_2);
        assertFalse("Sine switches are not constant outside the bounds.", sf.constantOutsideBounds());
        assertTrue("Sine switches are valid outside the bounds.", sf.validOutsideBounds());
        assertEquals("Default sine switch max-zero-derivative should return 1", 1, sf.getHighestOrderZeroDerivative());
        assertTrue("Default sine switch should be equal unity with symmetric inputs", sf.symmetricToUnity());
        
        for (double x = 0; x <= 1.0; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double ax = a*x;
            double sinOf = FastMath.sin(ax);
            double cosOf = FastMath.cos(ax);
            
            double valAt = sf.valueAt(x);
            double trueVal = sinOf * sinOf;
            assertEquals(String.format("Value of default sine switch at %8.4g should be %8.4g, was %8.4g", x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = sf.firstDerivative(x);
            trueVal = 2.0 * a * sinOf * cosOf;
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 100.0*ulp(trueVal);
            assertEquals(String.format("First derivative of default sine switch at %8.4g should be %8.4g, was %8.4g", x, trueVal, derivAt), trueVal, derivAt, delta);
            
            double d2 = sf.secondDerivative(x);
            trueVal = 2.0 * a * a * ((cosOf * cosOf) - (sinOf * sinOf));
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 100.0*ulp(trueVal);
            assertEquals(String.format("Second derivative of default sine switch at %8.4g should be %8.4g, was %8.4g", x, trueVal, d2), trueVal, d2, delta);
        }
        
        logger.info(" Testing manually-constructed default sine-squared switch.");
        
        sf = new SquaredTrigSwitch(piOverTwo, false);
        standardTest(sf, MID_TOLERANCE);
        a = piOverTwo;
        assertEquals("Default sine switch zero bound != 0.0", 0.0, sf.getZeroBound(), ULP_ZERO_2);
        assertEquals("Default sine switch one bound != 1.0", 1.0, sf.getOneBound(), ULP_ONE_2);
        assertFalse("Sine switches are not constant outside the bounds.", sf.constantOutsideBounds());
        assertTrue("Sine switches are valid outside the bounds.", sf.validOutsideBounds());
        assertEquals("Default sine switch max-zero-derivative should return 1", 1, sf.getHighestOrderZeroDerivative());
        assertTrue("Default sine switch should be equal unity with symmetric inputs", sf.symmetricToUnity());
        
        for (double x = 0; x <= 1.0; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double ax = a*x;
            double sinOf = FastMath.sin(ax);
            double cosOf = FastMath.cos(ax);
            
            double valAt = sf.valueAt(x);
            double trueVal = sinOf * sinOf;
            assertEquals(String.format("Value of default sine switch at %8.4g should be %8.4g, was %8.4g", x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = sf.firstDerivative(x);
            trueVal = 2.0 * a * sinOf * cosOf;
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 100.0*ulp(trueVal);
            assertEquals(String.format("First derivative of default sine switch at %8.4g should be %8.4g, was %8.4g", x, trueVal, derivAt), trueVal, derivAt, delta);
            
            double d2 = sf.secondDerivative(x);
            trueVal = 2.0 * a * a * ((cosOf * cosOf) - (sinOf * sinOf));
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 100.0*ulp(trueVal);
            assertEquals(String.format("Second derivative of default sine switch at %8.4g should be %8.4g, was %8.4g", x, trueVal, d2), trueVal, d2, delta);
        }
        
        logger.info(" Testing default cosine-squared switch.");
        
        sf = new SquaredTrigSwitch(true);
        standardTest(sf, MID_TOLERANCE);
        a = piOverTwo;
        assertEquals("Default cosine switch zero bound != 1.0", 1.0, sf.getZeroBound(), ULP_ONE_2);
        assertEquals("Default cosine switch one bound != 0.0", 0.0, sf.getOneBound(), ULP_ZERO_2);
        assertFalse("Cosine switches are not constant outside the bounds.", sf.constantOutsideBounds());
        assertTrue("Cosine switches are valid outside the bounds.", sf.validOutsideBounds());
        assertEquals("Default cosine switch max-zero-derivative should return 1", 1, sf.getHighestOrderZeroDerivative());
        assertTrue("Default cosine switch should be equal unity with symmetric inputs", sf.symmetricToUnity());
        
        for (double x = 0; x <= 1.0; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double ax = a*x;
            double sinOf = FastMath.sin(ax);
            double cosOf = FastMath.cos(ax);
            
            double valAt = sf.valueAt(x);
            double trueVal = cosOf * cosOf;
            assertEquals(String.format("Value of default cosine switch at %8.4g should be %8.4g, was %8.4g", x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = sf.firstDerivative(x);
            trueVal = -2.0 * a * sinOf * cosOf;
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 100.0*ulp(trueVal);
            assertEquals(String.format("First derivative of default cosine switch at %8.4g should be %8.4g, was %8.4g", x, trueVal, derivAt), trueVal, derivAt, delta);
            
            double d2 = sf.secondDerivative(x);
            trueVal = 2.0 * a * a * ((sinOf * sinOf) - (cosOf * cosOf));
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 100.0*ulp(trueVal);
            assertEquals(String.format("Second derivative of default cosine switch at %8.4g should be %8.4g, was %8.4g", x, trueVal, d2), trueVal, d2, delta);
        }
        
        logger.info(" Testing manually constructed default cosine-squared switch.");
        
        sf = new SquaredTrigSwitch(piOverTwo, true);
        standardTest(sf, MID_TOLERANCE);
        a = piOverTwo;
        assertEquals("Default cosine switch zero bound != 1.0", 1.0, sf.getZeroBound(), ULP_ONE_2);
        assertEquals("Default cosine switch one bound != 0.0", 0.0, sf.getOneBound(), ULP_ZERO_2);
        assertFalse("Cosine switches are not constant outside the bounds.", sf.constantOutsideBounds());
        assertTrue("Cosine switches are valid outside the bounds.", sf.validOutsideBounds());
        assertEquals("Default cosine switch max-zero-derivative should return 1", 1, sf.getHighestOrderZeroDerivative());
        assertTrue("Default cosine switch should be equal unity with symmetric inputs", sf.symmetricToUnity());
        
        for (double x = 0; x <= 1.0; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double ax = a*x;
            double sinOf = FastMath.sin(ax);
            double cosOf = FastMath.cos(ax);
            
            double valAt = sf.valueAt(x);
            double trueVal = cosOf * cosOf;
            assertEquals(String.format("Value of default cosine switch at %8.4g should be %8.4g, was %8.4g", x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = sf.firstDerivative(x);
            trueVal = -2.0 * a * sinOf * cosOf;
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 100.0*ulp(trueVal);
            assertEquals(String.format("First derivative of default cosine switch at %8.4g should be %8.4g, was %8.4g", x, trueVal, derivAt), trueVal, derivAt, delta);
            
            double d2 = sf.secondDerivative(x);
            trueVal = 2.0 * a * a * ((sinOf * sinOf) - (cosOf * cosOf));
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 100.0*ulp(trueVal);
            assertEquals(String.format("Second derivative of default cosine switch at %8.4g should be %8.4g, was %8.4g", x, trueVal, d2), trueVal, d2, delta);
        }
        
        logger.info(" Testing sine-squared switch with unadjusted (pi/2) bounds.");
        
        a = 1.0;
        sf = new SquaredTrigSwitch(a, false);
        standardTest(sf, MID_TOLERANCE);
        assertEquals(String.format("Sine switch %s zero bound != 0.0", sf), 0.0, sf.getZeroBound(), ULP_ZERO_2);
        assertEquals(String.format("Sine switch %s one bound != 1.0", sf), piOverTwo, sf.getOneBound(), ULP_ONE_2);
        assertFalse("Sine switches are not constant outside the bounds.", sf.constantOutsideBounds());
        assertTrue("Sine switches are valid outside the bounds.", sf.validOutsideBounds());
        assertEquals("Sine switch max-zero-derivative should return 1", 1, sf.getHighestOrderZeroDerivative());
        assertTrue("Sine switch should be equal unity with symmetric inputs", sf.symmetricToUnity());
        
        for (double x = 0; x <= 1.0; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double ax = a*x;
            double sinOf = FastMath.sin(ax);
            double cosOf = FastMath.cos(ax);
            
            double valAt = sf.valueAt(x);
            double trueVal = sinOf * sinOf;
            assertEquals(String.format("Value of sine switch %s at %8.4g should be %8.4g, was %8.4g", sf, x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = sf.firstDerivative(x);
            trueVal = 2.0 * a * sinOf * cosOf;
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 100.0*ulp(trueVal);
            assertEquals(String.format("First derivative of sine switch %s at %8.4g should be %8.4g, was %8.4g", sf, x, trueVal, derivAt), trueVal, derivAt, delta);
            
            double d2 = sf.secondDerivative(x);
            trueVal = 2.0 * a * a * ((cosOf * cosOf) - (sinOf * sinOf));
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 100.0*ulp(trueVal);
            assertEquals(String.format("Second derivative of sine switch %s at %8.4g should be %8.4g, was %8.4g", sf, x, trueVal, d2), trueVal, d2, delta);
        }
        
        logger.info(" Testing  sine-squared switch with doubled (2.0) bounds.");
        
        a = 0.5 * piOverTwo;
        sf = new SquaredTrigSwitch(a, false);
        standardTest(sf, MID_TOLERANCE);
        assertEquals(String.format("Sine switch %s zero bound != 0.0", sf), 0.0, sf.getZeroBound(), ULP_ZERO_2);
        assertEquals(String.format("Sine switch %s one bound != 1.0", sf), 2.0, sf.getOneBound(), 2.0 * ulp(2.0));
        assertFalse("Sine switches are not constant outside the bounds.", sf.constantOutsideBounds());
        assertTrue("Sine switches are valid outside the bounds.", sf.validOutsideBounds());
        assertEquals("Sine switch max-zero-derivative should return 1", 1, sf.getHighestOrderZeroDerivative());
        assertTrue("Sine switch should be equal unity with symmetric inputs", sf.symmetricToUnity());
        
        for (double x = 0; x <= 1.0; x += 0.01) {
            double delta = 10.0 * ulp(x);
            double ax = a*x;
            double sinOf = FastMath.sin(ax);
            double cosOf = FastMath.cos(ax);
            
            double valAt = sf.valueAt(x);
            double trueVal = sinOf * sinOf;
            assertEquals(String.format("Value of sine switch %s at %8.4g should be %8.4g, was %8.4g", sf, x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = sf.firstDerivative(x);
            trueVal = 2.0 * a * sinOf * cosOf;
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 100.0*ulp(trueVal);
            assertEquals(String.format("First derivative of sine switch %s at %8.4g should be %8.4g, was %8.4g", sf, x, trueVal, derivAt), trueVal, derivAt, delta);
            
            double d2 = sf.secondDerivative(x);
            trueVal = 2.0 * a * a * ((cosOf * cosOf) - (sinOf * sinOf));
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 100.0*ulp(trueVal);
            assertEquals(String.format("Second derivative of sine switch %s at %8.4g should be %8.4g, was %8.4g", sf, x, trueVal, d2), trueVal, d2, delta);
        }
        
        logger.info(" Testing cosine-squared switch with unadjusted (pi/2) bounds.");
        
        a = 1.0;
        sf = new SquaredTrigSwitch(a, true);
        standardTest(sf, MID_TOLERANCE);
        assertEquals(String.format("Cosine switch %s zero bound != 1.0", sf), piOverTwo, sf.getZeroBound(), 2.0*ulp(piOverTwo));
        assertEquals(String.format("Cosine switch %s one bound != 0.0", sf), 0.0, sf.getOneBound(), ULP_ZERO_2);
        assertFalse("Cosine switches are not constant outside the bounds.", sf.constantOutsideBounds());
        assertTrue("Cosine switches are valid outside the bounds.", sf.validOutsideBounds());
        assertEquals("Cosine switch max-zero-derivative should return 1", 1, sf.getHighestOrderZeroDerivative());
        assertTrue("Cosine switch should be equal unity with symmetric inputs", sf.symmetricToUnity());
        
        for (double x = 0; x <= 1.0; x += 0.01) {
            double delta = 50.0 * ulp(x);
            double ax = a*x;
            double sinOf = FastMath.sin(ax);
            double cosOf = FastMath.cos(ax);
            
            double valAt = sf.valueAt(x);
            double trueVal = cosOf * cosOf;
            assertEquals(String.format("Value of cosine switch %s at %8.4g should be %8.4g, was %8.4g", sf, x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = sf.firstDerivative(x);
            trueVal = -2.0 * a * sinOf * cosOf;
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 200.0*ulp(trueVal);
            assertEquals(String.format("First derivative of cosine switch %s at %8.4g should be %8.4g, was %8.4g", sf, x, trueVal, derivAt), trueVal, derivAt, delta);
            
            double d2 = sf.secondDerivative(x);
            trueVal = 2.0 * a * a * ((sinOf * sinOf) - (cosOf * cosOf));
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 200.0*ulp(trueVal);
            assertEquals(String.format("Second derivative of cosine switch %s at %8.4g should be %8.4g, was %8.4g", sf, x, trueVal, d2), trueVal, d2, delta);
        }
        
        logger.info(" Testing cosine-squared switch with doubled (2.0) bounds..");
        
        a = 0.5 * piOverTwo;
        sf = new SquaredTrigSwitch(a, true);
        standardTest(sf, MID_TOLERANCE);
        assertEquals(String.format("Cosine switch %s zero bound != 2.0", sf), 2.0, sf.getZeroBound(), 2.0*ulp(2.0));
        assertEquals(String.format("Cosine switch %s one bound != 0.0", sf), 0.0, sf.getOneBound(), ULP_ZERO_2);
        assertFalse("Cosine switches are not constant outside the bounds.", sf.constantOutsideBounds());
        assertTrue("Cosine switches are valid outside the bounds.", sf.validOutsideBounds());
        assertEquals("Cosine switch max-zero-derivative should return 1", 1, sf.getHighestOrderZeroDerivative());
        assertTrue("Cosine switch should be equal unity with symmetric inputs", sf.symmetricToUnity());
        
        for (double x = 0; x <= 1.0; x += 0.01) {
            double delta = 50.0 * ulp(x);
            double ax = a*x;
            double sinOf = FastMath.sin(ax);
            double cosOf = FastMath.cos(ax);
            
            double valAt = sf.valueAt(x);
            double trueVal = cosOf * cosOf;
            assertEquals(String.format("Value of cosine switch %s at %8.4g should be %8.4g, was %8.4g", sf, x, trueVal, valAt), trueVal, valAt, delta);
            
            double derivAt = sf.firstDerivative(x);
            trueVal = -2.0 * a * sinOf * cosOf;
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 200.0*ulp(trueVal);
            assertEquals(String.format("First derivative of cosine switch %s at %8.4g should be %8.4g, was %8.4g", sf, x, trueVal, derivAt), trueVal, derivAt, delta);
            
            double d2 = sf.secondDerivative(x);
            trueVal = 2.0 * a * a * ((sinOf * sinOf) - (cosOf * cosOf));
            delta = (trueVal < 1.0E-10) ? 1.0E-14 : 200.0*ulp(trueVal);
            assertEquals(String.format("Second derivative of cosine switch %s at %8.4g should be %8.4g, was %8.4g", sf, x, trueVal, d2), trueVal, d2, delta);
        }
    }
    
    /**
     * Standard set of tests that all implementations of 
     * UnivariateSwitchingFunction should pass; by default uses tight tolerances.
     * @param sf Switching function to test.
     */
    private void standardTest(UnivariateSwitchingFunction sf) {
        standardTest(sf, ULP_ONE_100);
    }
    
    /**
     * Standard set of tests that all implementations of 
     * UnivariateSwitchingFunction should pass. If looseTolerances is set, uses
     * a much looser tolerance for acceptance (1/1 million) instead of the 
     * default tolerances, approximately 100*ulp(0) and 100*ulp(1).
     * 
     * @param sf Switching function to test.
     * @param looseTolerances Use looser tolerances for test acceptance
     */
    private void standardTest(UnivariateSwitchingFunction sf, double tolerance) {
        double oneBound = sf.getOneBound();
        double zeroBound = sf.getZeroBound();
        double increment = ((oneBound - zeroBound) * 0.01);
        
        double minBound = 0.0 - tolerance;
        double maxBound = 1.0 + tolerance;
        
        for (int i = 0; i < 101; i++) {
            double pastLB = i * increment;
            // Slightly convoluted logic to avoid any round-off errors at ub.
            double x = (i == 100) ? oneBound : zeroBound + pastLB;
            double val = sf.valueAt(x);
            assertTrue(String.format("Switching function %s value at %8.4g was %8.4g, not in the range 0-1 inclusive", sf.toString(), x, val), val >= minBound && val <= maxBound);
            
            if (sf.symmetricToUnity()) {
                double symmX = oneBound - pastLB;
                double symmVal = sf.valueAt(symmX);
                assertEquals(String.format("Switching function %s should be "
                        + "symmetrical; values %7.4f and %7.4f at %7.4f and %7.4f "
                        + "do not sum to unity", sf.toString(), val, symmVal, x, symmX), 
                        1.0, (val + symmVal), tolerance);
            }
        }
        
        boolean validOutside = sf.validOutsideBounds();
        boolean constantOutside = sf.constantOutsideBounds();
        double valAtUB = sf.valueAt(sf.getOneBound());
        double valAtLB = sf.valueAt(sf.getZeroBound());
        if (Math.abs(valAtLB) < tolerance) {
            assertEquals(String.format("Switching function %s value at zero bound %8.4g was not 0.0 or 1.0, was %8.4g", sf.toString(), zeroBound, valAtLB), 0.0, valAtLB, tolerance);
            assertEquals(String.format("Switching function %s value at one bound %8.4g was not 0.0 or 1.0, was %8.4g", sf.toString(), oneBound, valAtUB), 1.0, valAtUB, tolerance);
        } else {
            assertEquals(String.format("Switching function %s value at zero bound %8.4g was not 0.0 or 1.0, was %8.4g", sf.toString(), zeroBound, valAtLB), 1.0, valAtLB, tolerance);
            assertEquals(String.format("Switching function %s value at one bound %8.4g was not 0.0 or 1.0, was %8.4g", sf.toString(), oneBound, valAtUB), 0.0, valAtUB, tolerance);
            logger.info(String.format(" Value of switching function %s at zero bound was 1.0, not 0.0; switching functions usually start at 0", sf));
        }
        if (validOutside || constantOutside) {
            for (int i = 1; i < 251; i++) {
                double pastBounds = i * increment;
                double x = zeroBound - pastBounds;
                double val = sf.valueAt(x);
                assertTrue(String.format("Switching function %s value at %8.4g (outside lb-ub) was %8.4g, not in the range 0-1 inclusive", sf.toString(), x, val), val > minBound && val < maxBound);
                if (constantOutside) {
                    assertEquals(String.format("Switching function %s value at %8.4g was %8.4g, did not match zero bound value %8.4g", sf.toString(), x, val, valAtLB), valAtLB, val, tolerance);
                }
                
                x = oneBound + pastBounds;
                val = sf.valueAt(x);
                assertTrue(String.format("Switching function %s value at %8.4g (outside lb-ub) was %8.4g, not in the range 0-1 inclusive", sf.toString(), x, val), val >= 0.0 && val <= 1.0);
                if (constantOutside) {
                    assertEquals(String.format("Switching function %s value at %8.4g was %8.4g, did not match one bound value %8.4g", sf.toString(), x, val, valAtUB), valAtUB, val, tolerance);
                }
            }
        }
        
        int maxZeroOrder = sf.getHighestOrderZeroDerivative();
        if (maxZeroOrder >= 1) {
            double deriv = sf.firstDerivative(zeroBound);
            assertEquals(String.format("Switching function %s first derivative at lb %8.4g was nonzero value %8.4g", sf.toString(), zeroBound, deriv), 0.0, deriv, tolerance);
            deriv = sf.firstDerivative(oneBound);
            assertEquals(String.format("Switching function %s first derivative at ub %8.4g was nonzero value %8.4g", sf.toString(), oneBound, deriv), 0.0, deriv, tolerance);
            if (maxZeroOrder >= 2) {
                deriv = sf.secondDerivative(zeroBound);
                assertEquals(String.format("Switching function %s second derivative at lb %8.4g was nonzero value %8.4g", sf.toString(), zeroBound, deriv), 0.0, deriv, tolerance);
                deriv = sf.secondDerivative(oneBound);
                assertEquals(String.format("Switching function %s second derivative at ub %8.4g was nonzero value %8.4g", sf.toString(), oneBound, deriv), 0.0, deriv, tolerance);
                for (int i = 3; i <= maxZeroOrder; i++) {
                    deriv = sf.nthDerivative(zeroBound, i);
                    assertEquals(String.format("Switching function %s %d-order derivative at lb %8.4g was nonzero value %8.4g", sf.toString(), i, zeroBound, deriv), 0.0, deriv, tolerance);
                    deriv = sf.nthDerivative(oneBound, i);
                    assertEquals(String.format("Switching function %s %d-order derivative at ub %8.4g was nonzero value %8.4g", sf.toString(), i, oneBound, deriv), 0.0, deriv, tolerance);
                }
            }
        }
        
        int maxOrderAtZero = sf.highestOrderZeroDerivativeAtZeroBound();
        if (maxOrderAtZero != maxZeroOrder && maxOrderAtZero > 0) {
            double deriv = sf.firstDerivative(zeroBound);
            assertEquals(String.format("Switching function %s first derivative at lb %8.4g was nonzero value %8.4g", sf.toString(), zeroBound, deriv), 0.0, deriv, tolerance);
            if (maxOrderAtZero >= 2) {
                deriv = sf.secondDerivative(zeroBound);
                assertEquals(String.format("Switching function %s second derivative at lb %8.4g was nonzero value %8.4g", sf.toString(), zeroBound, deriv), 0.0, deriv, tolerance);
                for (int i = 3; i <= maxZeroOrder; i++) {
                    deriv = sf.nthDerivative(zeroBound, i);
                    assertEquals(String.format("Switching function %s %d-order derivative at lb %8.4g was nonzero value %8.4g", sf.toString(), i, zeroBound, deriv), 0.0, deriv, tolerance);
                }
            }
        }
    }
}
