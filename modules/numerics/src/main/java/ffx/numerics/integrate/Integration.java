/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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

import java.text.DecimalFormat;

/**
 * This program integrates using three methods: the trapezoidal method,
 * Simpson's Three Point Integration, and Boole's Five Point Integration
 *
 * @author Claire O'Connell
 */
public class Integration {

    private final static double[] x = new double[202];

    static {
        x[0] = 0;
        for (int i = 1; i < 201; i++) {
            x[i] = .0025 + .005 * (i - 1);
        }
        x[201] = 1;
    }

    /**
     * <p>main.</p>
     *
     * @param args an array of {@link java.lang.String} objects.
     */
    public static void main(String[] args) {
        double testAnswer;
        double testTrap, testTrapRight, avgTrap, avgTrapError;
        double testSimp, testSimpRight, avgSimp, avgSimpError;
        double testBoole, testBooleRight, avgBoole, avgBooleError;
        double testRect, testRectRight, avgRect, avgRectError;
        DecimalFormat decimalFormat = new DecimalFormat("#.00");

        testAnswer = 2.98393938659796;
        System.out.print("The test case answer is " + testAnswer + "\n\n");

        testRect = rectangularMethodLeft(generateTestData_v1());
        //System.out.print("Testing rectangular left " + testRect + "\n");
        testRectRight = rectangularMethodRight(generateTestData_v1());
        //System.out.print("Testing rectangular right " + testRectRight + "\n");
        avgRect = (testRect + testRectRight) / 2.0;
        //System.out.print("Average rectangular " + avgRect + "\n");
        avgRectError = Math.abs(testAnswer - avgRect) / testAnswer * 100;
        //System.out.print("Average rectangular error " + decimalFormat.format(avgRectError) + "%\n");

        testTrap = trapInputLeft(generateTestData_v1());
        //System.out.print("Testing trapezoid left " + testTrap + "\n");
        testTrapRight = trapInputRight(generateTestData_v1());
        //System.out.print("Testing trapezoid right " + testTrapRight + "\n");
        avgTrap = (testTrap + testTrapRight) / 2.0;
        //System.out.print("Average trap " + avgTrap + "\n");
        avgTrapError = Math.abs(avgTrap - testAnswer) / testAnswer * 100;
        //System.out.print("Avg Trapezoidal error  " + decimalFormat.format(avgTrapError) + "%\n");

        testSimp = SimpsonsLeft(generateTestData_v1()) + HalfBinComposite(generateTestData_v1(), 1, "left");
        //System.out.print("Testing Simpsons left " + testSimp + "\n");
        testSimpRight = SimpsonsRight(generateTestData_v1()) + HalfBinComposite(generateTestData_v1(), 1, "right");
        //System.out.print("Testing Simpsons right " + testSimpRight + "\n");
        avgSimp = (testSimp + testSimpRight) / 2.0;
        //System.out.print("Average Simpsons " + avgSimp + "\n");
        avgSimpError = Math.abs(testAnswer - avgSimp) / testAnswer * 100;
        //System.out.print("Average Simpsons error " + decimalFormat.format(avgSimpError) + "%\n");

        testBoole = BooleLeft(generateTestData_v1()) + HalfBinComposite(generateTestData_v1(), 2, "left");
        //System.out.print("Testing Boole left " + testBoole + "\n");
        testBooleRight = BooleRight(generateTestData_v1()) + HalfBinComposite(generateTestData_v1(), 2, "right");
        //System.out.print("Testing Boole Right " + testBooleRight + "\n");
        avgBoole = (testBoole + testBooleRight) / 2.0;
        //System.out.print("Average Boole " + avgBoole + "\n");
        avgBooleError = Math.abs(testAnswer - avgBoole) / testAnswer * 100;
        //System.out.print("Average Boole error " + decimalFormat.format(avgBooleError) + "%\n");

        //Average integrals
        System.out.print("Average integrals \n\n");
        System.out.print("Average rectangular " + avgRect + "\n");
        System.out.print("Average trap " + avgTrap + "\n");
        System.out.print("Average Simpsons " + avgSimp + "\n");
        System.out.print("Average Boole " + avgBoole + "\n");

        //Average integral errors
        System.out.print("\nAverage integral error \n\n");
        System.out.print("Average rectangular error " + decimalFormat.format(avgRectError) + "%\n");
        System.out.print("Avg Trapezoidal error  " + decimalFormat.format(avgTrapError) + "%\n");
        System.out.print("Average Simpsons error " + decimalFormat.format(avgSimpError) + "%\n");
        System.out.print("Average Boole error " + decimalFormat.format(avgBooleError) + "%\n");
    }

    /**
     * <p>averageIntegral.</p>
     *
     * @param leftInt  a double.
     * @param rightInt a double.
     * @return a double.
     */
    public static double averageIntegral(double leftInt, double rightInt) {
        double avgInt = 0;

        avgInt = (leftInt + rightInt) / 2.0;
        //System.out.print("Average integral " + avgInt + "\n");
        return avgInt;
    }

    /**
     * <p>generateTestData_v1.</p>
     *
     * @return an array of {@link double} objects.
     */
    public static double[] generateTestData_v1() {
        double y[] = new double[202];

        for (int i = 0; i < 202; i++) {
            y[i] = 10 * Math.sin(6 * x[i]) - 7 * Math.cos(5 * x[i]) + 11 * Math.sin(8 * x[i]);
        }

        return y;
    }

    /**
     * <p>HalfBinComposite.</p>
     *
     * @param inputData an array of {@link double} objects.
     * @param mode      a int.
     * @param side      a {@link java.lang.String} object.
     * @return a double.
     */
    public static double HalfBinComposite(double[] inputData, int mode, String side) {
        double halfBinComposite = 0, lowerHalfBin, upperHalfBin, upperTrapArea;
        double lowerTrapArea, upperSimpson, lowerSimpson;
        int n;

        n = inputData.length;

        //Split by side first, then integration type
        if (side == "left") {
            //using trapezoidal integral for lower half bin
            lowerHalfBin = (inputData[1] + inputData[0]) / 2.0 * (x[1] - x[0]);
            halfBinComposite += lowerHalfBin;
            //System.out.print("lower bin " + halfBinComposite + "\n");
            switch (mode) {
                //case 1 is the Simpson's method, uses trapezoid on right for bin left out of Simpson's
                case 1:
                    upperTrapArea = (inputData[n - 3] + inputData[n - 2]) / 2.0 * (x[n - 2] - x[n - 3]);
                    //System.out.print("upper trap area " + upperTrapArea + "\n");
                    halfBinComposite += upperTrapArea;
                    break;
                //case 2 is the Boole's method, uses Simpsons and trapezoidal integral on right to cover remaining bins
                case 2:
                    upperSimpson = (1.0 / 3.0) * (x[n - 4] - x[n - 5]) * (inputData[n - 5] + 4 * inputData[n - 4] + inputData[n - 3]);
                    halfBinComposite += upperSimpson;
                    //System.out.print("Upper Simpson " + upperSimpson + "\n");
                    upperTrapArea = (inputData[n - 3] + inputData[n - 2]) / 2.0 * (x[n - 2] - x[n - 3]);
                    //System.out.print("upper trap area " + upperTrapArea + "\n");
                    halfBinComposite += upperTrapArea;
                    break;
            }
        } else if (side == "right") {
            //upper half bin calculated with trapezoid
            upperHalfBin = (inputData[n - 1] + inputData[n - 2]) / 2.0 * (x[n - 1] - x[n - 2]);
            halfBinComposite += upperHalfBin;
            switch (mode) {
                //case 1 is the Simpson's method, uses trapezoid on left for bin left out of Simpson's
                case 1:
                    lowerTrapArea = (inputData[1] + inputData[2]) / 2.0 * (x[2] - x[1]);
                    halfBinComposite += lowerTrapArea;
                    break;
                //case 2 is the Boole's method, uses Simpsons and trapezoidal integral on left to cover remaining bins
                case 2:
                    lowerTrapArea = (inputData[1] + inputData[2]) / 2.0 * (x[2] - x[1]);
                    halfBinComposite += lowerTrapArea;
                    lowerSimpson = (1.0 / 3.0) * (x[3] - x[2]) * (inputData[2] + 4 * inputData[3] + inputData[4]);
                    halfBinComposite += lowerSimpson;
                    break;
            }
        }
        //System.out.print("Half bin composite " + halfBinComposite + "\n");
        return halfBinComposite;
    }

    /**
     * <p>trapInputLeft.</p>
     *
     * @param inputData an array of {@link double} objects.
     * @return a double.
     */
    public static double trapInputLeft(double[] inputData) {
        double trapIntegral = 0, sum, area, total;

        int n = 0;

        n = x.length;

        sum = 0;
        total = 0;
        for (int a = 0; a < n - 2; a++) {
            if (a > 0) {
                area = (inputData[a + 1] + inputData[a]) / (double) 2 * (x[a + 1] - x[a]);
                sum = area + total;
                total = sum;
            }
            if (a == n - 3) {
                trapIntegral = sum;
                //System.out.print("\nThe trapezoidal integral is " + trapIntegral + "\n");
            }
            if (a == 0) {
                area = (inputData[a + 1] + inputData[a]) / (double) 2 * (x[a + 1] - x[a]);
                total = area;
            }
        }

        return trapIntegral;
    }

    /**
     * <p>trapInputRight.</p>
     *
     * @param inputData an array of {@link double} objects.
     * @return a double.
     */
    public static double trapInputRight(double[] inputData) {
        double trapIntegral = 0, sum, area, total;
        int n = 0;

        n = x.length;

        sum = 0;
        total = 0;
        for (int a = 1; a < n - 1; a++) {
            if (a > 1) {
                area = (inputData[a + 1] + inputData[a]) / 2.0 * (x[a + 1] - x[a]);
                sum = area + total;
                total = sum;
            }
            if (a == n - 2) {
                trapIntegral = sum;
                //System.out.print("\nThe trapezoidal integral is " + trapIntegral + "\n");

            }
            if (a == 1) {
                area = (inputData[a + 1] + inputData[a]) / 2.0 * (x[a + 1] - x[a]);
                total = area;
            }
        }

        return trapIntegral;
    }

    /**
     * <p>SimpsonsLeft.</p>
     *
     * @param inputData an array of {@link double} objects.
     * @return a double.
     */
    public static double SimpsonsLeft(double[] inputData) {
        double normalSimpsons = 0, area, sum, total;
        int n;

        n = inputData.length;

        sum = 0;
        total = 0;
        for (int a = 1; a < n - 4; a += 2) {
            area = (1.0 / 3.0) * (x[a + 1] - x[a]) * (inputData[a] + 4 * inputData[a + 1] + inputData[a + 2]);
            normalSimpsons += area;
                /*if (a == n-5){
                    System.out.print("Last Simpson's left starts at n-5 \n");
                }
                */
            //With half bin, goes into n-5
        }

        return normalSimpsons;
    }

    /**
     * <p>SimpsonsRight.</p>
     *
     * @param inputData an array of {@link double} objects.
     * @return a double.
     */
    public static double SimpsonsRight(double[] inputData) {
        double normalSimpsons = 0, area;
        int n;

        n = inputData.length;

        for (int a = 2; a < n - 3; a += 2) { //extra trap on lower edge so right edge of rightmost bin aligns with the upper half bin
            area = (1.0 / 3.0) * (x[a + 1] - x[a]) * (inputData[a] + 4 * inputData[a + 1] + inputData[a + 2]);
            normalSimpsons += area;
                /*if (a == n-4){ 
                    System.out.print("Last Simpson's right starts at n-4\n");
                }
                */

        }

        return normalSimpsons;
    }

    /**
     * <p>BooleLeft.</p>
     *
     * @param inputData an array of {@link double} objects.
     * @return a double.
     */
    public static double BooleLeft(double[] inputData) {
        double normalBoole = 0, area;
        int n;

        n = inputData.length;

        for (int a = 1; a < n - 5; a += 4) {
            area = (2.0 / 45.0) * (x[a + 1] - x[a]) * (7 * inputData[a] + 32 * inputData[a + 1] + 12 * inputData[a + 2] + 32 * inputData[a + 3] + 7 * inputData[a + 4]);
            normalBoole += area;
            
            /*if (a == n - 9) { //interval not compatible with 201 bins because of the half sized bin on the end
                System.out.print("Last Boole's left starts at n-9 \n"); 
            }
            */
        }

        return normalBoole;
    }

    /**
     * <p>BooleRight.</p>
     *
     * @param inputData an array of {@link double} objects.
     * @return a double.
     */
    public static double BooleRight(double[] inputData) {
        double normalBoole = 0, area;
        int n;

        n = inputData.length;


        for (int a = 4; a < n - 5; a += 4) { //Simpsons and trapezoid + lower bin on left
            area = (2.0 / 45.0) * (x[a + 1] - x[a]) * (7 * inputData[a] + 32 * inputData[a + 1] + 12 * inputData[a + 2] + 32 * inputData[a + 3] + 7 * inputData[a + 4]);
            normalBoole += area;
            
            /*if (a == n - 6) { //extra Simpsons on left edge so right edge of rightmost bin aligns with the upper half bin
                System.out.print("Last Boole's right starts at n-6 \n"); 
            }
                    */
        }

        return normalBoole;
    }

    /**
     * <p>rectangularMethodLeft.</p>
     *
     * @param inputData an array of {@link double} objects.
     * @return a double.
     */
    public static double rectangularMethodLeft(double[] inputData) {
        double rectangularIntegral = 0, area;
        double[] y = new double[202];
        int n;

        n = inputData.length;

        y = generateTestData_v1();

        for (int a = 0; a < n - 2; a++) {
            area = (x[a + 1] - x[a]) * y[a];
            rectangularIntegral += area;
        }

        //System.out.print("The left rectangular method is " + rectangularIntegral + "\n");
        return rectangularIntegral;
    }

    /**
     * <p>rectangularMethodRight.</p>
     *
     * @param inputData an array of {@link double} objects.
     * @return a double.
     */
    public static double rectangularMethodRight(double[] inputData) {
        double rectangularIntegral = 0, area = 0;
        double[] y = new double[202];
        int n;

        n = inputData.length;

        y = generateTestData_v1();

        for (int a = 1; a < n - 1; a++) {
            area = (x[a + 1] - x[a]) * y[a];
            rectangularIntegral += area;
        }

        //System.out.print("The right rectangular method is " + rectangularIntegral + "\n");
        return rectangularIntegral;
    }
}
