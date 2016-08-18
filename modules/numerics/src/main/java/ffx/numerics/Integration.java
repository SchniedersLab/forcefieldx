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

import static java.lang.Math.cos;
import static java.lang.Math.sin;

/**
 * This program integrates using three methods: the trapezoidal method,
 * Simpson's Three Point Integration, and Boole's Five Point Integration
 * @author ceoconnell
 */

public class Integration {

    public static void main(String[] args){  
        System.out.print("Rough integral approximations and errors for y=10sin(6x)-7cos(5x)+11sin(8x)\n");
        trapIntegral(1,201,.25); 
        SimpsonIntegral(1,201,.25); 
        BooleIntegral(1,201,.25); 
        
        System.out.print("\nSmooth integral approximations and errors for y=10sin(6x)-7cos(5x)+11sin(8x)\n");
        trapIntegral(1,201,.1); 
        SimpsonIntegral(1,201,.1); 
        BooleIntegral(1,201,.1); 
        
        //trapIntegral(1,201,.05); 
        //SimpsonIntegral(1,201,.05); 
        //BooleIntegral(1,201,.05); 
    }
    
    public static double trapIntegral(int lowerBound, int upperBound, double intervalSize){
        int a, n, i;
        double arraySize = (upperBound-lowerBound)*(1.0/intervalSize)+1;
        int size = (int) Math.ceil(arraySize);
        double sum, total, area, trapInt = 0, trapError;
        double[] y = new double [size];
        double[] x = new double [size];
        
        //Make sure that the interval is compatible with the trapezoidal method
        //Array size must be integer
        if (size != arraySize){
            System.out.print("Interval not compatible with trapezoidal integration\n");
            return -0;
        }
        
        for (i=0; i<(size);i++){
            x[i]=lowerBound+(double)intervalSize*i;
            y[i]=10*sin(6*x[i])-7*cos(5*x[i])+11*sin(8*x[i]);
        }
        n = x.length;
        
        a=0;
        sum=0;
        total=0;
        for (a=0;a<n-1;a++){
            if (a > 0){
                area = (y[a+1]+y[a])/(double)2*(x[a+1]-x[a]); 
                sum = area + total;
                total = sum;   
            }
            if (a == n-2){
                trapInt = sum;
                System.out.print("\nThe trapezoidal integral is " + trapInt + "\n");
                trapError = -2.2783 - trapInt;
                System.out.print("Trapezoidal error = " + trapError + "\n");
                 
            }
            if (a == 0){
                area = (y[a+1]+y[a])/(double)2*(x[a+1]-x[a]);
                total = area;
                //System.out.print("Total " + total + "\n");
            }
        }
       return trapInt; 
}
   
    public static double SimpsonIntegral(int lowerBound, int upperBound, double intervalSize){
            int a, n, i;
            double arraySize = (upperBound-lowerBound)*(1.0/intervalSize)+1;
            int size = (int) Math.ceil(arraySize);
            double sum, total, area, SimpsonInt = 0, SimpsonError;
            double[] y = new double [size];
            double[] x = new double [size];

            //Make sure that the interval is compatible with the Simpson's method
            //Array size must be multiple of 2n+1
            if (size != arraySize){
                System.out.print("Interval not compatible with Simpson's Three Point Integration\n");
                return -0;
            }

            for (i=0; i<(size);i++){
                x[i]=lowerBound+intervalSize*i;
                y[i]=10*sin(6*x[i])-7*cos(5*x[i])+11*sin(8*x[i]);
            }
            n = x.length;

            a=0;
            sum=0;
            total=0;
            for (a=0;a<n-2;a+=2){
                if (a > 0){
                    area = (1.0/3.0)*(x[a+1]-x[a])*(y[a]+4*y[a+1]+y[a+2]); 
                    sum = area + total;
                    total = sum;    
                }
                if (a == n-3){
                    SimpsonInt = sum;
                    System.out.print("The Simpson's Three Point integral is " + SimpsonInt + "\n");
                    SimpsonError = -2.2783 - SimpsonInt;
                    System.out.print("Simpson's error = " + SimpsonError + "\n");
                }
                if (a == 0){
                    area = (1.0/3.0)*(x[a+1]-x[a])*(y[a]+4*y[a+1]+y[a+2]); 
                    total = area;
                }
            }

       return SimpsonInt;    
    }

    public static double BooleIntegral(int lowerBound, int upperBound, double intervalSize){
            int a, n, i;
            double arraySize = (upperBound-lowerBound)*(1.0/intervalSize)+1;
            int size = (int) Math.ceil(arraySize);
            double sum, total, area, BooleInt = 0, BooleError;
            double[] y = new double [size];
            double[] x = new double [size];

            //Make sure interval is compatible with the five point method
            //Array size must be multiple of 4n+1
            if (size != arraySize){
                System.out.print("Interval not compatible with Boole's Five Point Integration\n");
                return -0;
            }

            for (i=0; i<(size);i++){
                x[i]=lowerBound+intervalSize*i;
                y[i]=10*sin(6*x[i])-7*cos(5*x[i])+11*sin(8*x[i]);
            }
            n = x.length;

            a=0;
            sum=0;
            total=0;
            for (a=0;a<n-4;a+=4){
                if (a > 0){
                    area = (2.0/45.0)*(x[a+1]-x[a])*(7*y[a]+32*y[a+1]+12*y[a+2]+32*y[a+3]+7*y[a+4]); 
                    sum = area + total;
                    total = sum;  
                }
                if (a == n-5){
                    BooleInt = sum;
                    System.out.print("The Boole's Three Point integral is " + BooleInt + "\n");
                    BooleError = -2.2783 - BooleInt;
                    System.out.print("Boole's error = " + BooleError + "\n");
                }
                if (a == 0){
                    area = (2.0/45.0)*(x[a+1]-x[a])*(7*y[a]+32*y[a+1]+12*y[a+2]+32*y[a+3]+7*y[a+4]); 
                    total = area;
                }
            }

       return BooleInt;    
    }   
}
