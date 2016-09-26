/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.numerics;

import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static ffx.numerics.Integration.HalfBinComposite;
import static ffx.numerics.Integration.generateTestData_v1;
import static ffx.numerics.Integration.trapInputLeft;
import static ffx.numerics.Integration.SimpsonsLeft;
import static ffx.numerics.Integration.BooleLeft;
import static ffx.numerics.Integration.rectangularMethodLeft;
import static ffx.numerics.Integration.trapInputRight;
import static ffx.numerics.Integration.SimpsonsRight;
import static ffx.numerics.Integration.BooleRight;
import static ffx.numerics.Integration.rectangularMethodRight;

/**The IntegrationTest is a JUnit test for the Integration program that ensures
 * that the integrals are calculated correctly. This test is run using known 
 * integrals calculated with the equation y=10sin(6x)-7cos(5x)+11sin(8x).
 * 
 * @author ceoconnell
 */
public class IntegrationTest {
    
    /**
    * Create array with pointers to doubles that will contain known integrals.
    */
    private double[] knownIntegral;
    
    /*
    Set the delta value for the assertEquals comparison
     */
    private final double DELTA = 1e-8;

    /**
     * Initializes the array before testing.
     */
    @Before
    public void setUp() {
        // Instantiate the knownIntegral array.
        knownIntegral = new double[8];
        
        /*The answers are in the order of the trapezoidal integral first, the
        Simpson's second, Boole's third, and rectangular method last. The 
        integrals are calculated with the bounds 1 and 201 with an interval of 
        .1. The first four are using left hand integrals and the second four
        use right hand integrals.
        */
        knownIntegral[0] = 2.9684353512887753;
        knownIntegral[1] = 2.9687126459508564;
        knownIntegral[2] = 2.968712622691869;
        knownIntegral[3] = 2.936215172510247;
        
        knownIntegral[4] = 3.0006927996084642;
        knownIntegral[5] = 3.000977174918476;
        knownIntegral[6] = 3.000977149598709;
        knownIntegral[7] = 2.968898509509485;
    }

    /**
    * Compares the calculated integrals with the known values.
    */
    @Test
    public void integrationTest() {

        /**Calculate the integrals using the left hand trapezoidal, Simpson's, 
         * and Boole's methods using data generated with the bounds 1 and 201 
         * with an interval of .1. The second four are the right handed integrals
         * in the same order.
        */ 
        double[] calculatedIntegral = new double[8];
       
        calculatedIntegral[0] = trapInputLeft(generateTestData_v1());
        calculatedIntegral[1] = SimpsonsLeft(generateTestData_v1())+HalfBinComposite(generateTestData_v1(),1,"left");
        calculatedIntegral[2] = BooleLeft(generateTestData_v1())+HalfBinComposite(generateTestData_v1(),2,"left");
        calculatedIntegral[3] = rectangularMethodLeft(generateTestData_v1());
        
        calculatedIntegral[4] = trapInputRight(generateTestData_v1());
        calculatedIntegral[5] = SimpsonsRight(generateTestData_v1())+HalfBinComposite(generateTestData_v1(),1,"right");
        calculatedIntegral[6] = BooleRight(generateTestData_v1())+HalfBinComposite(generateTestData_v1(),2,"right");
        calculatedIntegral[7] = rectangularMethodRight(generateTestData_v1());
        
        // Assert that the known integrals and calculated integrals are the same.
        for (int i=0;i<8;i++){
        assertEquals(knownIntegral[i],calculatedIntegral[i], DELTA);
        }
    }
}