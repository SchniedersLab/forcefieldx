//******************************************************************************
//
// File:    MSHistogramSeq.java
// Package: edu.rit.smp.fractal
// Unit:    Class edu.rit.smp.fractal.MSHistogramSeq
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

package edu.rit.smp.fractal;

//import edu.rit.pj.Comm;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

/**
 * Class MSHistogramSeq is a sequential program that calculates a histogram of
 * the Mandelbrot Set.
 * <P>
 * Usage: java edu.rit.smp.fractal.MSHistogramSeq <I>width</I> <I>height</I>
 * <I>xcenter</I> <I>ycenter</I> <I>resolution</I> <I>maxiter</I> <I>outfile</I>
 * <BR><I>width</I> = Image width (pixels)
 * <BR><I>height</I> = Image height (pixels)
 * <BR><I>xcenter</I> = X coordinate of center point
 * <BR><I>ycenter</I> = Y coordinate of center point
 * <BR><I>resolution</I> = Pixels per unit
 * <BR><I>maxiter</I> = Maximum number of iterations
 * <BR><I>outfile</I> = Output file name
 * <P>
 * The program calculates an iteration count <I>i</I> for each pixel in the same
 * way as the {@linkplain MandelbrotSetSeq} program, with 0 &lt;= <I>i</I> &lt;=
 * <I>maxiter</I>. The program prints into the output file a histogram of the
 * <I>i</I> values; that is, for each value of <I>i</I>, how many times that
 * value occurred in the image.
 * <P>
 * The computation is performed sequentially in a single processor. The program
 * measures the computation's running time, including the time to print the
 * output. This establishes a benchmark for measuring the computation's running
 * time on a parallel processor.
 *
 * @author  Alan Kaminsky
 * @version 02-Feb-2008
 */
public class MSHistogramSeq
	{

// Prevent construction.

	private MSHistogramSeq()
		{
		}

// Program shared variables.

	// Command line arguments.
	static int width;
	static int height;
	static double xcenter;
	static double ycenter;
	static double resolution;
	static int maxiter;
	static File outfile;

	// Initial pixel offsets from center.
	static int xoffset;
	static int yoffset;

	// Histogram (array of counters indexed by pixel value).
	static int[] histogram;

// Main program.

	/**
	 * Mandelbrot Set main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		//Comm.init (args);

		// Start timing.
		long t1 = System.currentTimeMillis();

		// Validate command line arguments.
		if (args.length != 7) usage();
		width = Integer.parseInt (args[0]);
		height = Integer.parseInt (args[1]);
		xcenter = Double.parseDouble (args[2]);
		ycenter = Double.parseDouble (args[3]);
		resolution = Double.parseDouble (args[4]);
		maxiter = Integer.parseInt (args[5]);
		outfile = new File (args[6]);

		// Initial pixel offsets from center.
		xoffset = -(width - 1) / 2;
		yoffset = (height - 1) / 2;

		// Create histogram.
		histogram = new int [maxiter + 1];

		long t2 = System.currentTimeMillis();

		// Compute all rows and columns.
		for (int r = 0; r < height; ++ r)
			{
			double y = ycenter + (yoffset - r) / resolution;

			for (int c = 0; c < width; ++ c)
				{
				double x = xcenter + (xoffset + c) / resolution;

				// Iterate until convergence.
				int i = 0;
				double aold = 0.0;
				double bold = 0.0;
				double a = 0.0;
				double b = 0.0;
				double zmagsqr = 0.0;
				while (i < maxiter && zmagsqr <= 4.0)
					{
					++ i;
					a = aold*aold - bold*bold + x;
					b = 2.0*aold*bold + y;
					zmagsqr = a*a + b*b;
					aold = a;
					bold = b;
					}

				// Increment histogram counter for pixel value.
				++ histogram[i];
				}
			}

		long t3 = System.currentTimeMillis();

		// Print histogram.
		PrintWriter out =
			new PrintWriter
				(new BufferedWriter
					(new FileWriter (outfile)));
		for (int i = 0; i <= maxiter; ++ i)
			{
			out.print (i);
			out.print ('\t');
			out.print (histogram[i]);
			out.println();
			}
		out.close();

		// Stop timing.
		long t4 = System.currentTimeMillis();
		System.out.println ((t2-t1) + " msec pre");
		System.out.println ((t3-t2) + " msec calc");
		System.out.println ((t4-t3) + " msec post");
		System.out.println ((t4-t1) + " msec total");
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.smp.fractal.MSHistogramSeq <width> <height> <xcenter> <ycenter> <resolution> <maxiter> <outfile>");
		System.err.println ("<width> = Image width (pixels)");
		System.err.println ("<height> = Image height (pixels)");
		System.err.println ("<xcenter> = X coordinate of center point");
		System.err.println ("<ycenter> = Y coordinate of center point");
		System.err.println ("<resolution> = Pixels per unit");
		System.err.println ("<maxiter> = Maximum number of iterations");
		System.err.println ("<outfile> = Output file name");
		System.exit (1);
		}

	}
