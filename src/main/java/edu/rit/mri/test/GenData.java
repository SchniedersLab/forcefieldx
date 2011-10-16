//******************************************************************************
//
// File:    GenData.java
// Package: edu.rit.mri.test
// Unit:    Class edu.rit.mri.test.GenData
//
// This Java source file is copyright (C) 2008 by Alan Kaminsky. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// ark@cs.rit.edu.
//
// This Java source file is part of the Parallel Java Library ("PJ"). PJ is free
// software; you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// PJ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************

package edu.rit.mri.test;

import edu.rit.numeric.NormalPrng;

import edu.rit.util.Random;

/**
 * Class GenData is a main program that generates a spin signal data set.
 * <P>
 * Usage: java edu.rit.mri.test.GenData <I>seed</I> <I>stdev</I> <I>tmax</I>
 * <I>nt</I> <I>&rho;</I><SUB>1</SUB> <I>x</I><SUB>1</SUB> [
 * <I>&rho;</I><SUB>2</SUB> <I>x</I><SUB>2</SUB> . . . ]
 * <P>
 * The program:
 * <OL TYPE=1>
 * <LI>
 * Initializes a pseudorandom number generator with <I>seed</I>.
 * <LI>
 * Generates <I>nt</I> values of <I>t</I>, unevenly spaced, up to <I>tmax</I>.
 * <LI>
 * For each value of <I>t</I>, computes the spin signal
 * <CENTER>
 * <I>S</I>(<I>t</I>)&emsp;=&emsp;&Sigma;<SUB><I>i</I></SUB>&emsp;<I>&rho;</I><SUB><I>i</I></SUB> [1 &minus; 2 exp(&minus;<I>x</I><SUB><I>i</I></SUB> <I>t</I>)]
 * </CENTER>
 * where <I>&rho;</I><SUB><I>i</I></SUB> is the <I>i</I>-th spin density and
 * <I>x</I><SUB><I>i</I></SUB> is the <I>i</I>-th spin relaxation rate.
 * <LI>
 * For each spin signal value, adds random Gaussian noise with mean 0 and
 * standard deviation <I>stdev</I>.
 * <LI>
 * Prints each value of <I>t</I> and the corresponding <I>S</I>(<I>t</I>) with
 * noise on one line of the standard output.
 * </OL>
 *
 * @author  Alan Kaminsky
 * @version 10-Jun-2008
 */
public class GenData
	{

// Prevent construction.

	private GenData()
		{
		}

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Parse command line arguments.
		if (args.length < 6 || (args.length % 2) != 0) usage();
		long seed = Long.parseLong (args[0]);
		double stdev = Double.parseDouble (args[1]);
		double tmax = Double.parseDouble (args[2]);
		int nt = Integer.parseInt (args[3]);
		int ns = (args.length - 4)/2;
		double[] rho = new double [ns];
		double[] x = new double [ns];
		for (int i = 0; i < ns; ++ i)
			{
			rho[i] = Double.parseDouble (args[2*i+4]);
			x[i] = Double.parseDouble (args[2*i+5]);
			}

		// Set up Gaussian PRNG.
		NormalPrng prng =
			new NormalPrng (Random.getInstance (seed), 0.0, stdev);

		// Generate and print data.
		for (int i = 1; i <= nt; ++ i)
			{
			//double t = Math.pow (tmax, ((double) i)/((double) nt));
			double t = Math.pow (((double) i)/((double) nt), 2.4) * tmax;
			double s = 0.0;
			for (int j = 0; j < ns; ++ j)
				{
				s += rho[j]*(1.0 - 2.0*Math.exp(-x[j]*t));
				}
			s += prng.next();
			System.out.print (t);
			System.out.print ('\t');
			System.out.print (s);
			System.out.println();
			}
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mri.test.GenData <seed> <stdev> <tmax> <nt> <rho1> <x1> [<rho2> <x2> ...]");
		System.err.println ("<seed> = Random seed");
		System.err.println ("<stdev> = Measurement error standard deviation");
		System.err.println ("<tmax> = Largest t value");
		System.err.println ("<nt> = Number of t values");
		System.err.println ("<rho1> = Spin density");
		System.err.println ("<x1> = Spin relaxation rate");
		System.exit (1);
		}

	}
