/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.potential.bonded;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Residue.ResidueType;
import ffx.potential.bonded.AminoAcidUtils;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.parsers.PDBFilter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.util.FastMath;
import java.util.logging.Logger;

/**
 *
 * @author mrtollefson
 */
public class SturmMethod {
    
    private static final Logger logger = Logger.getLogger(SturmMethod.class.getName());
    
    final int PRINT_LEVEL = 0;
    final int MAX_ORDER = 16;
    final int MAXPOW = 32;
    final double SMALL_ENOUGH = 1.0e-18;
    double RELERROR;
    int MAXIT, MAX_ITER_SECANT;
    double[] roots;
    
    /*Class that holds polynomial information.
    Poly objects have one position that holds the order of a polynomial in
          poly.ord
    Poly objects have up to 16 positions that hold the coefficients of a 
          16th degree polynomial. The coefficients are located 
          in poly.coef[position] */
    private class poly
    {
        int ord;
        double coef[];
        
        private poly() 
        {
            coef = new double[MAX_ORDER+1];
        }
    }

    //sets termination criteria for polynomial solver
    void initialize_sturm(double[][] tol_secant, int[][] max_iter_sturm, int[][] max_iter_secant) {
        RELERROR = tol_secant[0][0];
        MAXIT = max_iter_sturm[0][0];
        MAX_ITER_SECANT = max_iter_secant[0][0];
    }

    
    void solve_sturm(int[] p_order, int[] n_root, double[] poly_coeffs, double[] roots) {
        poly[] sseq = new poly[MAX_ORDER*2];
        double min, max;
        int order, i, j, nroots, nchanges, np;
        int[] atmin=new int[1];
        int[] atmax=new int[1];
        this.roots = roots;    
        /*for (i=0;i<16;i++)
        {
            System.out.println("roots: "+roots[i]);
        } */
        
        for (i=0; i<MAX_ORDER*2; i++)
        {
            sseq[i]=new poly();
        }
        order = p_order[0];
        //System.out.println("p_order[0]: "+p_order[0]+"\n");
        
        for (i = order; i >= 0; i--)
        {
            sseq[0].coef[i] = poly_coeffs[i];
            //System.out.println("sseq: "+sseq[0].coef[i]+"\n");
        }

       /* for (j=0; j<=32; j++)
        {
            for (i = 0; i<=16; i++)
            {
                System.out.println("first sseq[1]: "+sseq[j].coef[i]);
            }    
            System.out.println("\n\n\n");
        }*/
        if (PRINT_LEVEL > 0) 
        {
            StringBuilder string = new StringBuilder();
            for (i = order; i >= 0; i--) 
            {
                string.append(String.format("coefficients in Sturm solver\n"));
                string.append(String.format("%d %f\n", i, sseq[0].coef[i]));
            }
            logger.info(string.toString());
        }
        
        //System.out.println("SSEQ: "+Arrays.deepToString(sseq)+"\n\n\n\n");
        //build the sturm sequence
        //System.out.println("Order: "+order+"\n");
        /*int counter;
        int counter1;
        for (counter = 0; counter <32 ; counter++)
        {
            for (counter1 = 0; counter1 <16 ; counter1 ++)
            {
                System.out.println("sseq: "+ sseq[counter].coef[counter1]);
            }
        } */
        np = buildsturm(order, sseq);
        
        /*int counter;
        int counter1;
        for (counter = 0; counter <32 ; counter++)
        {
            for (counter1 = 0; counter1 <16 ; counter1 ++)
            {
                System.out.println("sseq: "+ sseq[counter].coef[counter1]);
            }
        } */
        
        //System.out.println("NP: "+np);
        if (PRINT_LEVEL > 0) 
        {   
            StringBuilder string1 = new StringBuilder();
            string1.append(String.format("Sturm sequence for:\n"));
            for (i = order; i >= 0; i--)
            {
                string1.append(String.format("%f ", sseq[0].coef[i]));
                string1.append("\n");
            }    
            for (i = 0; i <= np; i++) 
            {   
                for (j = sseq[i].ord; j >= 0; j--)
                {
                    string1.append(String.format("%f ", sseq[i].coef[j]));
                    string1.append("\n");
                }
            }
            logger.info(string1.toString());
	}

        //get the number of real roots
        nroots = numroots(np, sseq, atmin, atmax);
        //System.out.println("NROOTS: "+nroots+"\n");
                
        if (nroots == 0)
        {
            n_root[0] = nroots;
        }
        
        if (PRINT_LEVEL > 0)
        {
            logger.info(String.format("Number of real roots: %d\n", nroots));
        }
        
        //calculate the bracket that the roots live in
        min = -1.0;
        
        //System.out.println("np: "+np);
    //int counter;
    //int counter1;
    //for (counter = 0; counter <32 ; counter++)
    //{
    //    for (counter1 = 0; counter1 <16 ; counter1 ++)
    //    {
    //        System.out.println("sseq: "+ sseq[counter].coef[counter1]);
    //    }
    //}
        
        nchanges = numchanges(np, sseq, min);
        
        for (i = 0; nchanges != atmin[0] && i != MAXPOW; i++)
        { 
            min *= 10.0;
            nchanges = numchanges(np, sseq, min);
        }
        
        if (nchanges != atmin[0]) 
        {
            logger.info(String.format("solve: unable to bracket all negative roots\n"));
            atmin[0] = nchanges;
        }
        
        max = 1.0;
        nchanges = numchanges(np, sseq, max);
        
        for (i = 0; nchanges != atmax[0] && i != MAXPOW; i++) 
        { 
          max *= 10.0;
          nchanges = numchanges(np, sseq, max);
        }
        if (nchanges != atmax[0]) 
        {
          logger.info(String.format("solve: unable to bracket all positive roots\n"));
          atmax[0] = nchanges;
        }
        
        //perform the bisection
        nroots = atmin[0] - atmax[0];
        
        //ERROR!!!
        //for(i=0;i<16;i++)
        //{
        //    System.out.println("first_roots: " + roots[i]);
        //}

        sbisect(np, sseq, min, max, atmin[0], atmax[0], this.roots);
        
        //for(i=0;i<16;i++)
        //{
        //    System.out.println("roots!!!: " + this.roots[i]);
        //}
        
        n_root[0] = nroots;
        
        //write out the roots
        if (PRINT_LEVEL > 0) 
        {
            if (nroots == 1) 
            {
                logger.info(String.format("\n1 distinct real root at x = %f\n", this.roots[0]));
            } 
            else 
            {
                StringBuilder string2 = new StringBuilder();
                string2.append(String.format("\n%d distinct real roots for x: \n", nroots));
                for (i = 0; i != nroots; i++)
                  {
                    string2.append(String.format("%f\n", this.roots[i]));
                  }
                logger.info(string2.toString());
            }
        }
    }
    
    double hyper_tan(double a, double x) {
        double exp_x1, exp_x2, ax;
        
        ax = a*x;
        
        if (ax > 100.0)
        {
          return(1.0);
        }
        else if (ax < -100.0)
        {
          return(-1.0);
        }
        else     
        {
          exp_x1 = Math.exp(ax);
          exp_x2 = Math.exp(-ax);
          return (exp_x1 - exp_x2)/(exp_x1 + exp_x2);
        }
    }
    
    
    //calculates the modulus of u(x)/v(x) leaving it in r, it returns 0 if r(x)
    //      is constant. This function assumes the leading coefficient of v is
    //      is 1 or -1.
    static boolean modp(poly u, poly v, poly r) {
        
        //Modp was originally a static int and returned r.ord as an int. The
        //      buildsturm function requires that a boolean is returned
        //      from the modp function, so modp is set to return a boolean.
        //      This new boolean return could be problematic elsewhere in the
        //      code if modp is used to return a value. This should be checked.
        int k;
        int j;
	double[] nr; 
        double end;
        double[] uc;
        int i;
        
        nr = r.coef;
        end = u.ord;
        uc = u.coef;
        
        //System.out.println("uc: "+Arrays.toString(u.coef));
        for (i=0;i<=end;i++)
        {
            nr[i]=uc[i];
           // System.out.println("NR[I]:"+nr[i]);
        }
        //System.out.println("V: "+v.coef[v.ord]+"\n");
        if(v.coef[v.ord]<0.0)
        {
            for(k = u.ord-v.ord-1; k>=0; k-=2)
            {
                r.coef[k]=-r.coef[k];
                //System.out.println("r.coef[k]: "+r.coef[k]+"\n");
            }
            
            for(k = u.ord-v.ord; k>=0; k--)
            {
                for(j=v.ord+k-1; j>=k; j--)
                {
                    r.coef[j]=-r.coef[j] - r.coef[v.ord+k] * v.coef[j-k];
                }
            }
        }
        else
        {
            for (k=u.ord-v.ord; k>=0; k--)
            {
                for(j=v.ord+k-1; j>= k; j--)
                {
                    r.coef[j] -= r.coef[v.ord+k]*v.coef[j-k];
                    //System.out.println("r.coef[j]: "+r.coef[j]+"\n");
                    //System.out.println("K: "+k+"\n");
                    //System.out.println("J: "+j+"\n");
                }
            }
        }
        
        k=v.ord-1;
        //System.out.println("K: "+k+"\n");
        while(k>=0 && FastMath.abs(r.coef[k]) < 1.0e-18)
        {
            r.coef[k]=0.0;
            k--;
        }
        
        if(k<0)
        {
            r.ord = 0;
        }
        else
        {
            r.ord = k;
        }
        
        //returns the boolean answer
        if(r.ord>0)
        {
            //System.out.println("TRUE\n\n");
            return (true);
        }
        
        else if(r.ord==0)
        {
            //System.out.println("False\n\n");
            return (false);
        }
        
        else
        {
            //System.out.println("FALSE\n\n");
            return (false);
        }
    }
    
    //build up a sturm sequence for a polynomial, and return the number of
    //      polynomials in the sequence
    int buildsturm(int ord, poly[] sseq) {
        int i;
        int j;
	double f;
        double[] fp;
        double[] fc;
 
        sseq[0].ord=ord;
        sseq[1].ord=ord-1;
        
        //print statements for debugging purposes
        //System.out.println("sseq[0].ord: "+sseq[0].ord+"\n");
        //System.out.println("sseq[1].ord: "+sseq[1].ord+"\n");
   
        // calculate the derivative and normalise the leading coefficient
        f= Math.abs(sseq[0].coef[ord]*ord);
        //System.out.println("F: "+f+"\n");
        fp = sseq[1].coef;
        //System.out.println("fp: "+Arrays.toString(fp)+"\n");
        //System.out.println("sseq[1].coef: "+Arrays.toString(sseq[1].coef));
   
	fc = Arrays.copyOfRange(sseq[0].coef, 1, sseq[0].coef.length);
        //System.out.println("fc: "+Arrays.toString(fc)+"\n");
        //System.out.println("sseq[0].coef: "+sseq[0].coef[2]);
        //System.out.println("fp: "+Arrays.toString(fp)+"\n");
        
        j = 0;
        for (i = 1; i <= ord; i++)
        {
            //print statements for debugging purposes
            //System.out.println("*fc++: "+fc[j]);
            //System.out.println("i: "+i);
            //System.out.println("f: "+f);
            //System.out.println("i/f: "+i/f);
            //System.out.println("fc++*258: "+fc[j]*258.084);
            
            fp[j]=fc[j]*i/f;
            
            //print statement for debugging purposes
            //System.out.println("*fp++: "+fp[j]+"\n\n\n");
            j++;
        }

        /*System.out.println("fp: "+Arrays.toString(fp)+"\n");
        for( i = 0; i<16; i++)
        {
            System.out.println("second sseq[1]: "+sseq[1].coef[i]);
        }*/
        
        //construct the rest of the Sturm sequence 
        //System.out.println("sseq[i]:" + Arrays.toString(sseq) + "\n");
        for(i=0;i<sseq[0].coef.length-2 && modp(sseq[i],sseq[i+1],sseq[i+2]);i++)
        {   
            //reverse the sign and normalise
            f = -Math.abs(sseq[i+2].coef[sseq[i+2].ord]);
            //System.out.println("F: "+f+"\n");
            
            for(j=(int)sseq[i+2].ord;j>=0;j--)
            {
                sseq[i+2].coef[j] /= f;
                //fp[(int)j]/=f;
                //System.out.println("sseq: "+sseq[i+2].coef[j]+"\n");
            }
        }

        sseq[i+2].coef[0]=-sseq[i+2].coef[0]; //reverse the sign
        //System.out.println("spi+2.coef[0]: "+sseq[i+2].coef[0]);
        
        //System.out.println("sseq.length-i: " + (sseq.length-i) + "\n");
        //System.out.println(sseq[0].ord);
        //System.out.println(sseq[i+2].ord);
        //System.out.println("subtract: "+(sseq[0].ord-sseq[i+2].ord)+"\n");
        
        return(sseq[0].ord-sseq[i+2].ord);
    }
    
    //return the number of distinct real roots of the polynomial described in sseq
    int numroots(int np, poly[] sseq, int[] atneg, int[] atpos) {
        int atposinf, atneginf;
	double	f, lf;
        int i;
 
        atposinf = atneginf = 0;
       
        //changes at positive infinity
        lf = sseq[0].coef[sseq[0].ord];
        
        //Print statements for debugging purposes
        //System.out.println("The value of np is:\n");
        //System.out.print(np);
        
        //ERROR2!!!!!
        //for(i=0;i<sseq[0].coef.length-2 && modp(sseq[i],sseq[i+1],sseq[i+2]);i++)
        //for (sp = sseq + 2; modp(sp - 2, sp - 1, sp); sp++)
        
        //int counter = 0;
        for(i=1;i<=np;i++)
        //for (s = sseq + 1; s <= sseq + np; s++)
        {
            f=sseq[i].coef[sseq[i].ord];
            //System.out.println("F: "+f);
            if (lf == 0.0 || lf * f < 0)
            {
                atposinf++;
            }
            
            lf=f;
            //System.out.println("lf: "+lf);
            //System.out.println("atp: "+atposinf);
            //counter++;
        }
        //System.out.println("Counter: "+counter+"\n");
        
        //changes at negative infinity
        if ((sseq[0].ord & 1)!=0)
        {
            lf = -sseq[0].coef[sseq[0].ord];
        }
        else
        {
            lf = sseq[0].coef[sseq[0].ord];
        }
        
        for(i=1;i<=np;i++)
        {
            if((sseq[i].ord & 1)!=0)
            {
                f=-sseq[i].coef[sseq[i].ord];
            }
            else
            {
                f=sseq[i].coef[sseq[i].ord];
            }
            if (lf == 0.0 || lf * f < 0)
            {
                atneginf++;
            }
            
            lf = f;
            //System.out.println("lf: "+lf);
            //System.out.println("atneginf: "+atneginf);
        }
        
        atneg[0] = atneginf;
	atpos[0] = atposinf;
        
        //System.out.println("atneginf: "+(atneginf-atposinf));
        return (atneginf - atposinf);
    }
    
    //return the number of sign changesin the Sturm sequence in sseq 
    //      at the value a.
    int numchanges(int np, poly[] sseq, double a) {
        int changes;
        double f;
        double lf;
        //poly[] s;
        
        changes = 0;
        //ERROR!!!
        //System.out.println("sseq.ord: "+sseq[0].ord);
        //System.out.println("a: "+a);
        lf = evalpoly(sseq[0].ord, sseq[0].coef, a);
        //System.out.println("Lf: "+lf);
        for (int i=1; i<=np; i++) 
        {
            f = evalpoly(sseq[i].ord, sseq[i].coef, a);
            if (lf == 0.0 || lf * f < 0)
            {
		changes++;
            }
            
            lf = f;
        }
        
        //System.out.println(changes);
        return(changes);
    }
    
    //uses a bisection based on the sturm sequence for the polynomial described
    //      in sseq to isolate intervals in which roots occur, the roots are 
    //      returned in the roots array in order of magnitude
    
    void sbisect(int np, poly[] sseq, double min, double max, int atmin, int atmax, double[] roots) {
        double	mid=0.;
        int n1 = 0;
        int n2 = 0;
        int its, atmid, nroot;
        int remainder = 0;
        int count = 0;
        
        nroot = atmin-atmax;
        
        if (nroot==1)
            {
                //first try a less expensive technique
                if (modrf(sseq[0].ord, sseq[0].coef, min, max, roots))
                {
                    remainder = this.roots.length - roots.length;
                    count = 0;
                    for(int q = remainder; q<this.roots.length;q++){
                        this.roots[q] = roots[count];
                        count++;
                    }
                        return;
                    
                }
                
                
                //When code reaches this point, the root must be evaluated
                //      using the Sturm sequence.
                for(its=0; its<MAXIT; its++)
                {
                    mid = (min+max)/2;
                    
                    atmid = numchanges(np, sseq, mid);
                    
                    if (Math.abs(mid)>RELERROR)
                    {
                        if (Math.abs((max-min)/mid)<RELERROR)
                        {
                            roots[0]=mid;
                            
                            remainder = this.roots.length - roots.length;
                            count = 0;
                            for (int q = remainder; q < this.roots.length; q++) {
                                this.roots[q] = roots[count];
                                count++;
                            }
                            return;
                            
                        }
                    }
                    else if (Math.abs(max-min)<RELERROR)
                    {
                        roots[0]=mid;
                        
                        
                        remainder = this.roots.length - roots.length;
                        count = 0;
                        for (int q = remainder; q < this.roots.length; q++) {
                            this.roots[q] = roots[count];
                            count++;
                        }
                        return;

                    }
                    
                    if((atmin-atmid)==0)
                    {
                        min = mid;
                    }
                    else
                    {
                        max = mid;
                    }
                }
                if (its==MAXIT)
                {
                    logger.info(String.format("sbisect: overflow min %f max %f diff %f nroot %d n1 %d n2 %d\n", min, max, max - min, nroot, n1, n2));
                    roots[0]=mid;
                }
                
                remainder = this.roots.length - roots.length;
                count = 0;
                for (int q = remainder; q < this.roots.length; q++) {
                    this.roots[q] = roots[count];
                    count++;
                }
                return;

            }
        
        //more than one root in the interval, must bisect...
        for (its = 0; its < MAXIT; its++)
        {
            mid = (min + max) / 2;
      
            atmid = numchanges(np, sseq, mid);

            n1 = atmin - atmid;
            n2 = atmid - atmax;
            //System.out.println("N1: " +n1+"\n");
            //System.out.println("N2: " +n2+"\n");

            if (n1 != 0 && n2 != 0) 
              {
                  sbisect(np, sseq, min, mid, atmin, atmid, roots);
                  sbisect(np, sseq, mid, max, atmid, atmax, Arrays.copyOfRange(roots, n1, roots.length));
                  break;
              }

            if (n1 == 0)
            {
                min = mid;
            }
            else
            {
                max = mid;
            }
        }

        if (its == MAXIT) 
          {
              for(n1 = atmax; n1 < atmin; n1++)
              {
                  roots[n1 - atmax] = mid;
                  //System.out.println("rootsat: "+roots[n1-atmax]);
              }
          }
        
        //may not need this~~~~
        remainder = this.roots.length - roots.length;
        count = 0;
        for (int q = remainder; q < this.roots.length; q++) {
            this.roots[q] = roots[count];
            count++;
        }
    }
    
    //evaluate polynomial defined in coef returning its value
    double evalpoly(int ord, double[] coef, double x) {
        double[] fp=new double[ord];
        double f;
        
        fp = coef;
	f = fp[ord];

        int i = ord;
        for (i--;i>=0;i--)
        {
            f = x*f + fp[i];
        }
	return(f);
    }
    
    //uses the modified regula-falsi method to evaluate the root in interval
    //      [a,b] of the polynomial described in coef. The root is returned in
    //      val. The routine returns zero if it can't converge.
    boolean modrf(int ord, double[] coef, double a, double b, double[] val) {
        int its;
        int ecoef;
        double fa, fb, x, fx, lfx;
        double[] fp, scoef;
        int i;
        
        scoef=coef;
        ecoef=ord;
        
        fa=fb=coef[ord];
        for(i=ord-1;i>=0;i--)
        {
         fa = a * fa + coef[i];
         fb = b * fb + coef[i];   
        }
        
        if (fa * fb > 0.0)
        {
            return(false);
        }
        
        lfx = fa;
        
        for (its = 0; its < MAX_ITER_SECANT; its++) 
        {
            x = (fb * a - fa * b) / (fb - fa);
            
            // constrain that x stays in the bounds
            if (x < a || x > b)
            {
                x = 0.5 * (a+b);
            }
            
            fx=coef[ord];
            
            for(i=ord-1;i>=0;i--)
            {
                fx=x*fx+coef[i];   
            }
            
            if (Math.abs(x) > RELERROR) 
            {
              if (Math.abs(fx / x) < RELERROR) 
                {
                  val[0] = x;
                  return(true);
                }
            }
             else if (Math.abs(fx) < RELERROR) 
            {
                val[0] = x;
                return(true);
            }
            if ((fa * fx) < 0) 
            {
              b = x;
              fb = fx;
              if ((lfx * fx) > 0)
              fa /= 2;
            } 
            else 
              {
                a = x;
                fa = fx;
                if ((lfx * fx) > 0)
                  fb /= 2;
              }

              lfx = fx;
        }
 
        return(false);   
    }
    
    File write_pdb_backbone(char[] pdb_name, String[] res_name, double[][] r_n, double[][] r_a, double[][] r_c, int stt_res, int end_res, char[] chain_n, char[] chain_a, char[] chain_c,MolecularAssembly molAss, int counter, boolean writeFile)
    {
        Polymer[] newChain = molAss.getChains();
        ArrayList<Atom> backBoneAtoms;
        double[] xyz_n = new double[3];
        double[] xyz_a = new double[3];
        double[] xyz_c = new double[3];
        //start position = start position of loop, same for end
        
        for(int i = stt_res + 1; i < end_res; i++){
            Residue newResidue = newChain[0].getResidue(i);
            backBoneAtoms = newResidue.getBackboneAtoms();
            for (Atom backBoneAtom : backBoneAtoms) {
                //System.out.println("ATOM:"+backBoneAtom.getAtomType().name);
                switch(backBoneAtom.getAtomType().name)
                {
                    case "C":
                        xyz_c[0] = r_c[i-stt_res][0]; 
                        xyz_c[1] = r_c[i-stt_res][1]; 
                        xyz_c[2] = r_c[i-stt_res][2]; 
                       //System.out.println("C coordinates: "+Arrays.toString(xyz_c)+"\n");
                        
                        backBoneAtom.moveTo(xyz_c);
                        
                        //System.out.println("C coordinates: "+Arrays.toString(backBoneAtom.getXYZ())+"\n");
                        break;
                    case "N":
                        xyz_n[0] = r_n[i-stt_res][0];
                        xyz_n[1] = r_n[i-stt_res][1];
                        xyz_n[2] = r_n[i-stt_res][2];
                        //molAss.getBackBoneAtoms().get(i).setXYZ(xyz_n);
                        //System.out.println("N coordinates: "+Arrays.toString(xyz_n)+"\n");
                        backBoneAtom.moveTo(xyz_n);
                        break;
                    case "CA":
                        xyz_a[0] = r_a[i-stt_res][0];
                        xyz_a[1] = r_a[i-stt_res][1];
                        xyz_a[2] = r_a[i-stt_res][2];
                        backBoneAtom.moveTo(xyz_a);
                        //System.out.println("A coordinates: "+Arrays.toString(xyz_a)+"\n");
                        //molAss.getBackBoneAtoms().get(i).setXYZ(xyz_a);
                        break;
                    case "HA":
                        newResidue.deleteAtom(backBoneAtom);
                        break;
                    case "H":
                        //System.out.println("ATOM:"+backBoneAtom.getAtomType().name);
                        newResidue.deleteAtom(backBoneAtom);
                        break;
                    default:
                        newResidue.deleteAtom(backBoneAtom);
                        break;
                }
            }
            
            ArrayList<Atom> sideChainAtoms = newResidue.getSideChainAtoms();
            for (Atom sideChainAtom:sideChainAtoms)
            {
                newResidue.deleteAtom(sideChainAtom);
            }
        }
 
        File file = molAss.getFile();
        
        for(int i = stt_res; i <= end_res; i++)
        {
           // if(molAss.getBackBoneAtoms().get(i).getAtomType().name.equals("C"))
          //  {
                double[] xyz = molAss.getBackBoneAtoms().get(i).getXYZ();
                //System.out.println(Arrays.toString(xyz)+"\n");
          //  }
        }
        
        String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
        if(!filename.contains("_loop")){
            filename = filename+"_loop";
        }
        
        File modifiedFile = new File(filename + ".pdb_"+counter);
        PDBFilter modFilter = new PDBFilter(modifiedFile,molAss,null,null);
        
        if(writeFile)
        {
            modFilter.writeFile(modifiedFile, true);
        }
        
        /*int res_no, ir, k=0;
        for(res_no=stt_res;res_no<=end_res;res_no++)
        {
            ir = res_no - stt_res;
            k++;
            System.out.print("ATOM\t"+k+" N\t"+res_name[ir]+" "+chain_n[ir]+"\t"+res_no+"\t");
            System.out.printf("%.3f %.3f %.3f \n",r_n[ir][0],r_n[ir][1],r_n[ir][2]);
            k++;
            System.out.printf("ATOM\t"+k+" CA\t"+res_name[ir]+" "+chain_a[ir]+"\t"+res_no+"\t");
            System.out.printf("%.3f %.3f %.3f \n",r_a[ir][0],r_a[ir][1],r_a[ir][2]);
            k++;
            System.out.printf("ATOM\t"+k+" C\t"+res_name[ir]+" "+chain_c[ir]+"\t"+res_no+"\t");
            System.out.printf("%.3f %.3f %.3f \n",r_c[ir][0],r_c[ir][1],r_c[ir][2]);
        }
        System.out.println("\n\n\n"); */
        
        return (modifiedFile);
    }
}
