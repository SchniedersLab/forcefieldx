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
package ffx.numerics;

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
        for (int i = 0; i < 202; i++) {
            x[i] = 0; //TODO 
        }
    }
    
    public static void main(String[] args) {
        double trapInt, BooleInt, SimpsonInt, rectInt, testAnswer;
        double trapError = 0, BooleError = 0, SimpsonError = 0, rectError = 0;
        double avgTrap, avgSimpson, avgBoole, avgRect;
        double trapErrorR, SimpsonErrorR, BooleErrorR, rectErrorR;
        double trapIntRight, SimpsonIntRight, BooleIntRight, rectIntRight;
        double avgTrapError, avgSimpsonError, avgBooleError, avgRectError;
        DecimalFormat decimalFormat = new DecimalFormat("#.00");

        System.out.print("The left handed integrals are: \n");
        trapInt = trapInputLeft(generateTestDataLeft_v1());
        SimpsonInt = SimpsonInputLeft(generateTestDataLeft_v1());
        BooleInt = BooleInputLeft(generateTestDataLeft_v1());
        rectInt = rectangularMethodLeft(generateTestDataLeft_v1());

        testAnswer = 2.98393938659796;
        trapError = Math.abs(trapInt - testAnswer) / testAnswer * 100;
        SimpsonError = Math.abs(SimpsonInt - testAnswer) / testAnswer * 100;
        BooleError = Math.abs(BooleInt - testAnswer) / testAnswer * 100;
        rectError = Math.abs(rectInt - testAnswer) / testAnswer * 100;

        System.out.print("\nLeft hand errors by integration method \n\n");
        System.out.print("Trapezoidal error  " + decimalFormat.format(trapError) + "%\n");
        System.out.print("Simpson error      " + decimalFormat.format(SimpsonError) + "%\n");
        System.out.print("Boole error        " + decimalFormat.format(BooleError) + "%\n");
        System.out.print("Rectangular error  " + decimalFormat.format(rectError) + "%\n");

        System.out.print("\nThe right handed integrals are: \n");
        trapIntRight = trapInputRight(generateTestDataRight_v1()); 
        SimpsonIntRight = SimpsonInputRight(generateTestDataRight_v1());
        BooleIntRight = BooleInputRight(generateTestDataRight_v1());
        rectIntRight = rectangularMethodRight(generateTestDataRight_v1());

        trapErrorR = Math.abs(trapIntRight - testAnswer) / testAnswer * 100;
        SimpsonErrorR = Math.abs(SimpsonIntRight - testAnswer) / testAnswer * 100;
        BooleErrorR = Math.abs(BooleIntRight - testAnswer) / testAnswer * 100;
        rectErrorR = Math.abs(rectIntRight - testAnswer) / testAnswer * 100;

        System.out.print("\nRight hand errors by integration method \n\n");
        System.out.print("Trapezoidal error  " + decimalFormat.format(trapErrorR) + "%\n");
        System.out.print("Simpson error      " + decimalFormat.format(SimpsonErrorR) + "%\n");
        System.out.print("Boole error        " + decimalFormat.format(BooleErrorR) + "%\n");
        System.out.print("Rectangular error  " + decimalFormat.format(rectErrorR) + "%\n");

        System.out.print("\nAverage Integrals\n\n");
        avgTrap = averageIntegral(trapInt, trapIntRight);
        avgSimpson = averageIntegral(SimpsonInt, SimpsonIntRight);
        avgBoole = averageIntegral(BooleInt, BooleIntRight);
        avgRect = averageIntegral(rectInt, rectIntRight);

        avgTrapError = Math.abs(avgTrap - testAnswer) / testAnswer * 100;
        avgSimpsonError = Math.abs(avgSimpson - testAnswer) / testAnswer * 100;
        avgBooleError = Math.abs(avgBoole - testAnswer) / testAnswer * 100;
        avgRectError = Math.abs(avgRect - testAnswer) / testAnswer * 100;

        System.out.print("\nAveraged errors by integration method \n\n");
        System.out.print("Trapezoidal error  " + decimalFormat.format(avgTrapError) + "%\n");
        System.out.print("Simpson error      " + decimalFormat.format(avgSimpsonError) + "%\n");
        System.out.print("Boole error        " + decimalFormat.format(avgBooleError) + "%\n");
        System.out.print("Rectangular error  " + decimalFormat.format(avgRectError) + "%\n");
    }

    public static double averageIntegral(double leftInt, double rightInt) {
        double avgInt = 0;

        avgInt = (leftInt + rightInt) / 2.0;
        System.out.print("Average integral " + avgInt + "\n");
        return avgInt;
    }

    public static double[] generateTestDataLeft_v1() {
        double y[] = new double[201];
        double x[] = new double[201];

        x[0] = 0;
        x[1] = .0025;
        for (int i = 2; i < 201; i++) {
            x[i] = .0025 + .005 * (i - 1);
        }

        for (int i = 0; i < 201; i++) {
            y[i] = 10 * Math.sin(6 * x[i]) - 7 * Math.cos(5 * x[i]) + 11 * Math.sin(8 * x[i]);
        }

        return y;
    }

    public static double[] generateTestDataRight_v1() {
        double y[] = new double[201];
        double x[] = new double[201];

        for (int i = 0; i < 200; i++) {
            x[i] = .0025 + .005 * i;
        }
        x[200] = 1.0;

        for (int i = 0; i < 201; i++) {
            y[i] = 10 * Math.sin(6 * x[i]) - 7 * Math.cos(5 * x[i]) + 11 * Math.sin(8 * x[i]);
        }

        return y;
    }

    public static double trapInputLeft(double[] inputData) {
        double trapIntegral = 0, sum, area, total;
        double[] x = new double[201];
        int n = 0;

        x[0] = 0;
        x[1] = .0025;
        //Initializes the normal sized bins 
        for (int i = 2; i < 201; i++) {
            x[i] = .0025 + .005 * (i - 1);
        }

        /*for (int j=0;j<201;j++){
         System.out.print("x[" + j + "]= " + x[j] + "\n");
         }
         */
        n = x.length;

        sum = 0;
        total = 0;
        for (int a = 0; a < n - 1; a++) {
            if (a > 0) {
                area = (inputData[a + 1] + inputData[a]) / (double) 2 * (x[a + 1] - x[a]);
                sum = area + total;
                total = sum;
            }
            if (a == n - 2) {
                trapIntegral = sum;
                System.out.print("\nThe trapezoidal integral is " + trapIntegral + "\n");
            }
            if (a == 0) {
                area = (inputData[a + 1] + inputData[a]) / (double) 2 * (x[a + 1] - x[a]);
                total = area;
            }
        }

        return trapIntegral;
    }

    public static double SimpsonInputLeft(double[] inputData) {
        double SimpsonIntegral = 0, SimpsonError, sum, area = 0, total;
        double lowerTrapArea, upperTrapArea;
        double[] x = new double[201];
        int n;

        x[0] = 0;
        x[1] = .0025;
        //Initializes the normal sized bins 
        for (int i = 2; i < 201; i++) {
            x[i] = .0025 + .005 * (i - 1);
        }

        /*for (int j=0;j<201;j++){
         System.out.print("x[" + j + "]= " + x[j] + "\n");
         }     
         */
        n = x.length;

        sum = 0;
        total = 0;
        for (int a = 1; a < n - 2; a += 2) {
            if (a > 0) {
                area = (1.0 / 3.0) * (x[a + 1] - x[a]) * (inputData[a] + 4 * inputData[a + 1] + inputData[a + 2]);
                sum = area + total;
                total = sum;
            }
            if (a == n - 4) {
                SimpsonIntegral = sum;
            }

        }

        lowerTrapArea = (inputData[0] + inputData[1]) / (double) 2 * (x[1] - x[0]);
        SimpsonIntegral += lowerTrapArea;
        upperTrapArea = (inputData[n - 2] + inputData[n - 1]) / (double) 2 * (x[n - 1] - x[n - 2]);
        SimpsonIntegral += upperTrapArea;

        System.out.print("The Simpson's Three Point integral with two trapezoidal integrals is " + SimpsonIntegral + "\n");

        return SimpsonIntegral;
    }

    public static double BooleInputLeft(double[] inputData) {
        double BooleIntegral = 0, sum, area = 0, total;
        double lowerTrapArea, upperTrapArea, BooleSimpsonIntegral;
        double upperSimpsonArea, upperTrapArea3 = 0;
        double[] x = new double[201];
        int n;

        x[0] = 0;
        x[1] = .0025;
        //Initializes the normal sized bins 
        for (int i = 2; i < 201; i++) {
            x[i] = .0025 + .005 * (i - 1);
        }
        /*
         for (int j=0;j<201;j++){
         System.out.print("x[" + j + "]= " + x[j] + "\n");
         }
         */
        n = x.length;

        sum = 0;
        total = 0;
        for (int a = 1; a < n - 4; a += 4) {
            if (a > 0) {
                area = (2.0 / 45.0) * (x[a + 1] - x[a]) * (7 * inputData[a] + 32 * inputData[a + 1] + 12 * inputData[a + 2] + 32 * inputData[a + 3] + 7 * inputData[a + 4]);
                sum = area + total;
                total = sum;
            }
            if (a == n - 8) { //interval not compatible with 201 bins because of the half sized bin on the end
                BooleIntegral = sum;
            }
        }

        BooleSimpsonIntegral = BooleIntegral;
        lowerTrapArea = (inputData[0] + inputData[1]) / (double) 2 * (x[1] - x[0]);
        BooleIntegral += lowerTrapArea;
        BooleSimpsonIntegral += lowerTrapArea;

        //Calculates the remaining upper bins with one Simpson's and one trapezoidal integral
        upperSimpsonArea = (1.0 / 3.0) * (x[n - 3] - x[n - 4]) * (inputData[n - 4] + 4 * inputData[n - 3] + inputData[n - 2]);
        upperTrapArea = (inputData[n - 2] + inputData[n - 1]) / (double) 2 * (x[n - 1] - x[n - 2]);
        BooleSimpsonIntegral += (upperSimpsonArea + upperTrapArea);
        System.out.print("The Boole's Five Point integral with one Simpson and two trapezoidal integrals is " + BooleSimpsonIntegral + "\n");

        //Calculates the remaining upper bins using three trapezoidal integrals
        for (int a = 0; a < 3; a++) {
            area = (inputData[a + 198] + inputData[a + 197]) / (double) 2 * (x[a + 198] - x[a + 197]);
            upperTrapArea3 += area;
        }
        BooleIntegral += upperTrapArea3;
        System.out.print("The Boole's Five Point integral with four trapezoidal integrals is " + BooleIntegral + "\n");

        return BooleSimpsonIntegral;
    }

    public static double rectangularMethodLeft(double[] inputData) {
        double rectangularIntegral = 0, sum = 0, area = 0;
        double[] x = new double[201];
        double[] y = new double[201];
        int n;
        x[0] = 0;
        x[1] = .0025;
        //Initializes the normal sized bins 
        for (int i = 2; i < 201; i++) {
            x[i] = .0025 + .005 * (i - 1);
        }

        y = generateTestDataLeft_v1();

        for (int a = 0; a < 200; a++) {
            area = (x[a + 1] - x[a]) * y[a];
            rectangularIntegral += area;
        }

        System.out.print("The old integration method is " + rectangularIntegral + "\n");
        return rectangularIntegral;
    }

    public static double trapInputRight(double[] inputData) {
        double trapIntegral = 0, sum, area, total;
        double[] x = new double[201];
        int n = 0;

        //Initializes the normal sized bins 
        for (int i = 0; i < 201; i++) {
            x[i] = .0025 + .005 * i;
        }
        x[200] = 1.0;
        /*for (int j=0;j<201;j++){
         System.out.print("x[" + j + "]= " + x[j] + "\n");
         }
         */
        n = x.length;

        sum = 0;
        total = 0;
        for (int a = 0; a < n - 1; a++) {
            if (a > 0) {
                area = (inputData[a + 1] + inputData[a]) / (double) 2 * (x[a + 1] - x[a]);
                sum = area + total;
                total = sum;
            }
            if (a == n - 2) {
                trapIntegral = sum;
                System.out.print("\nThe trapezoidal integral is " + trapIntegral + "\n");

            }
            if (a == 0) {
                area = (inputData[a + 1] + inputData[a]) / (double) 2 * (x[a + 1] - x[a]);
                total = area;
            }
        }

        return trapIntegral;
    }

    //SIMPSON AND BOOLE RIGHT HAND METHODS ARE CURRENTLY NOT FUNCTIONAL :)
    public static double SimpsonInputRight(double[] inputData) {
        double SimpsonIntegral = 0, sum, area = 0, total;
        double upperTrapArea = 0;
        double[] x = new double[201];
        int n;

        //Initializes the normal sized bins 
        for (int i = 0; i < 200; i++) {
            x[i] = .0025 + .005 * i;
        }
        x[200] = 1.0;

        /*for (int j=0;j<201;j++){
         System.out.print("x[" + j + "]= " + x[j] + "\n");
         }
         */
        n = x.length;

        sum = 0;
        total = 0;
        for (int a = 0; a < n - 4; a += 2) {
            if (a > 0) {
                area = (1.0 / 3.0) * (x[a + 1] - x[a]) * (inputData[a] + 4 * inputData[a + 1] + inputData[a + 2]);
                sum = area + total;
                total = sum;
            }
            if (a == n - 5) {
                SimpsonIntegral = sum;
            }

        }
        for (int a = 0; a < 2; a++) {
            area = (inputData[n - 2 + a] + inputData[n - 3 + a]) / (double) 2 * (x[n - 2 + a] - x[n - 3 + a]);
            upperTrapArea += area;
        }
        SimpsonIntegral += upperTrapArea;

        System.out.print("The Simpson's Three Point integral is " + SimpsonIntegral + "\n");
                //System.out.print("Total " + total + "\n");

        return SimpsonIntegral;
    }

    public static double BooleInputRight(double[] inputData) {
        double BooleIntegral = 0, BooleError, sum, area = 0, total;
        double upperTrapArea2 = 0, BooleSimpsonIntegral, BooleTrapIntegral;
        double upperSimpsonArea, upperTrapArea4 = 0;
        double[] x = new double[201];
        int n;

        //Initializes the normal sized bins 
        for (int i = 0; i < 201; i++) {
            x[i] = .0025 + .005 * i;
        }
        x[200] = 1.0;
        /*for (int j=0;j<201;j++){
         System.out.print("x[" + j + "]= " + x[j] + "\n");
         }
         */
        n = x.length;

        sum = 0;
        total = 0;
        for (int a = 0; a < n - 6; a += 4) {
            if (a > 0) {
                area = (2.0 / 45.0) * (x[a + 1] - x[a]) * (7 * inputData[a] + 32 * inputData[a + 1] + 12 * inputData[a + 2] + 32 * inputData[a + 3] + 7 * inputData[a + 4]);
                sum = area + total;
                total = sum;
            }
            if (a == n - 9) {//interval not compatible with 201 bins
                BooleIntegral = sum;
                    //System.out.print("The Boole's Five Point integral is " + BooleIntegral + "\n");

            }
        }

        BooleSimpsonIntegral = BooleIntegral;
        BooleTrapIntegral = BooleIntegral;

        //Calculates the remaining upper bins with one Simpson's and two trapezoidal integral
        upperSimpsonArea = (1.0 / 3.0) * (x[n - 4] - x[n - 5]) * (inputData[n - 5] + 4 * inputData[n - 4] + inputData[n - 3]);
        for (int a = 0; a < 2; a++) {
            area = (inputData[n - 2] + inputData[n - 1]) / (double) 2 * (x[n - 1] - x[n - 2]);
            upperTrapArea2 += area;
        }

        BooleSimpsonIntegral += (upperSimpsonArea + upperTrapArea2);
        System.out.print("The Boole's Five Point integral with one Simpson and two trapezoidal integrals is " + BooleSimpsonIntegral + "\n");

        //Calculates the remaining upper bins using three trapezoidal integrals
        for (int a = 0; a < 4; a++) {
            area = (inputData[a + n - 4] + inputData[a + n - 5]) / (double) 2 * (x[a + n - 4] - x[a + n - 5]);
            upperTrapArea4 += area;
        }
        BooleTrapIntegral += upperTrapArea4;
        System.out.print("The Boole's Five Point integral with four trapezoidal integrals is " + BooleTrapIntegral + "\n");

        return BooleTrapIntegral;
    }

    public static double rectangularMethodRight(double[] inputData) {
        double rectangularIntegral = 0, sum = 0, area = 0;
        double[] x = new double[201];
        double[] y = new double[201];
        int n;

        //Initializes the normal sized bins 
        for (int i = 0; i < 201; i++) {
            x[i] = .0025 + .005 * i;
        }
        x[200] = 1.0;

        y = generateTestDataRight_v1();

        for (int a = 0; a < 200; a++) {
            area = (x[a + 1] - x[a]) * y[a];
            rectangularIntegral += area;
        }

        System.out.print("The old integration method is " + rectangularIntegral + "\n");
        return rectangularIntegral;
    }

}
