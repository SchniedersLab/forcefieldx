//******************************************************************************
//
// File:    MandelbrotSetSmp2.java
// Package: edu.rit.smp.fractal
// Unit:    Class edu.rit.smp.fractal.MandelbrotSetSmp2
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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

import edu.rit.color.HSB;

import edu.rit.image.PJGColorImage;
import edu.rit.image.PJGImage;

import edu.rit.pj.BarrierAction;
import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;

import java.util.ArrayList;

/**
 * Class MandelbrotSetSmp2 is an SMP parallel program that calculates the
 * Mandelbrot Set.
 * <P>
 * Usage: java -Dpj.nt=<I>K</I> edu.rit.smp.fractal.MandelbrotSetSmp2
 * <I>width</I> <I>height</I> <I>xcenter</I> <I>ycenter</I> <I>resolution</I>
 * <I>maxiter</I> <I>gamma</I> <I>filename</I>
 * <BR><I>K</I> = Number of processors
 * <BR><I>width</I> = Image width (pixels)
 * <BR><I>height</I> = Image height (pixels)
 * <BR><I>xcenter</I> = X coordinate of center point
 * <BR><I>ycenter</I> = Y coordinate of center point
 * <BR><I>resolution</I> = Pixels per unit
 * <BR><I>maxiter</I> = Maximum number of iterations
 * <BR><I>gamma</I> = Used to calculate pixel hues
 * <BR><I>filename</I> = PJG image file name
 * <P>
 * The program considers a rectangular region of the complex plane centered at
 * (<I>xcenter,ycenter</I>) of <I>width</I> pixels by <I>height</I> pixels,
 * where the distance between adjacent pixels is 1/<I>resolution</I>. The
 * program takes each pixel's location as a complex number <I>c</I> and performs
 * the following iteration:
 * <P>
 * <I>z</I><SUB>0</SUB> = 0
 * <BR><I>z</I><SUB><I>i</I>+1</SUB> = <I>z</I><SUB><I>i</I></SUB><SUP>2</SUP> + <I>c</I>
 * <P>
 * until <I>z</I><SUB><I>i</I></SUB>'s magnitude becomes greater than or equal
 * to 2, or <I>i</I> reaches a limit of <I>maxiter</I>. The complex numbers
 * <I>c</I> where <I>i</I> reaches a limit of <I>maxiter</I> are considered to
 * be in the Mandelbrot Set. (Actually, a number is in the Mandelbrot Set only
 * if the iteration would continue forever without <I>z</I><SUB><I>i</I></SUB>
 * becoming infinite; the foregoing is just an approximation.) The program
 * creates an image with the pixels corresponding to the complex numbers
 * <I>c</I> and the pixels' colors corresponding to the value of <I>i</I>
 * achieved by the iteration. Following the traditional practice, points in the
 * Mandelbrot set are black, and the other points are brightly colored in a
 * range of colors depending on <I>i</I>. The exact hue of each pixel is
 * (<I>i</I>/<I>maxiter</I>)<SUP><I>gamma</I></SUP>. The image is stored in a
 * Parallel Java Graphics (PJG) file specified on the command line.
 * <P>
 * The computation is performed in parallel in multiple processors. The program
 * measures the computation's running time, including the time to write the
 * image file. The program also measures the running time for each parallel
 * team thread and for each chunk of indexes computed by the parallel for loop.
 *
 * @author  Alan Kaminsky
 * @version 11-Nov-2007
 */
public class MandelbrotSetSmp2
	{

// Prevent construction.

	private MandelbrotSetSmp2()
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
	static double gamma;
	static File filename;

	// Initial pixel offsets from center.
	static int xoffset;
	static int yoffset;

	// Image matrix.
	static int[][] matrix;
	static PJGColorImage image;

	// Table of hues.
	static int[] huetable;

	// Starting and ending times for each thread.
	static long[] thr_t1;
	static long[] thr_t2;

	// Starting and ending times for each chunk of indexes in each thread.
	static ArrayList<Long>[] chunk_t1;
	static ArrayList<Long>[] chunk_t2;

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
		if (args.length != 8) usage();
		width = Integer.parseInt (args[0]);
		height = Integer.parseInt (args[1]);
		xcenter = Double.parseDouble (args[2]);
		ycenter = Double.parseDouble (args[3]);
		resolution = Double.parseDouble (args[4]);
		maxiter = Integer.parseInt (args[5]);
		gamma = Double.parseDouble (args[6]);
		filename = new File (args[7]);

		// Initial pixel offsets from center.
		xoffset = -(width - 1) / 2;
		yoffset = (height - 1) / 2;

		// Create image matrix to store results.
		matrix = new int [height] [width];
		image = new PJGColorImage (height, width, matrix);

		// Create table of hues for different iteration counts.
		huetable = new int [maxiter+1];
		for (int i = 0; i < maxiter; ++ i)
			{
			huetable[i] = HSB.pack
				(/*hue*/ (float) Math.pow (((double)i)/((double)maxiter),gamma),
				 /*sat*/ 1.0f,
				 /*bri*/ 1.0f);
			}
		huetable[maxiter] = HSB.pack (1.0f, 1.0f, 0.0f);

		// Starting and ending times for each thread and each index chunk.
		int K = ParallelTeam.getDefaultThreadCount();
		thr_t1 = new long [K];
		thr_t2 = new long [K];
		chunk_t1 = new ArrayList [K];
		chunk_t2 = new ArrayList [K];
		for (int i = 0; i < K; ++ i)
			{
			chunk_t1[i] = new ArrayList<Long>();
			chunk_t2[i] = new ArrayList<Long>();
			}

		long t2 = System.currentTimeMillis();

		// Compute all rows and columns.
		new ParallelTeam().execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				final int thr = getThreadIndex();
				final ArrayList<Long> chunk_t1_thr = chunk_t1[thr];
				final ArrayList<Long> chunk_t2_thr = chunk_t2[thr];
				thr_t1[thr] = System.currentTimeMillis();

				execute (0, height-1, new IntegerForLoop()
					{
					public void run (int first, int last)
						{
						chunk_t1_thr.add (System.currentTimeMillis());
						for (int r = first; r <= last; ++ r)
							{
							int[] matrix_r = matrix[r];
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

								// Record number of iterations for pixel.
								matrix_r[c] = huetable[i];
								}
							}
						chunk_t2_thr.add (System.currentTimeMillis());
						}
					},
				BarrierAction.NO_WAIT);

				thr_t2[thr] = System.currentTimeMillis();
				}
			});

		long t3 = System.currentTimeMillis();

		// Write image to PJG file.
		PJGImage.Writer writer =
			image.prepareToWrite
				(new BufferedOutputStream
					(new FileOutputStream (filename)));
		writer.write();
		writer.close();

		// Stop timing.
		long t4 = System.currentTimeMillis();
		System.out.println ((t2-t1) + " msec pre finish");
		for (int i = 0; i < K; ++ i)
			{
			System.out.println
				("\t" + (thr_t1[i]-t1) + " msec thread " + i + " start");
			for (int j = 0; j < chunk_t1[i].size(); ++ j)
				{
				System.out.println
					("\t\t" + (chunk_t1[i].get(j)-t1) + " msec chunk " + j +
					 " start");
				System.out.println
					("\t\t" + (chunk_t2[i].get(j)-t1) + " msec chunk " + j +
					 " finish");
				}
			System.out.println
				("\t" + (thr_t2[i]-t1) + " msec thread " + i + " finish");
			}
		System.out.println ((t3-t1) + " msec calc finish");
		System.out.println ((t4-t1) + " msec post finish");
		System.out.println ((t4-t1) + " msec total");
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.nt=<K> edu.rit.smp.fractal.MandelbrotSetSmp2 <width> <height> <xcenter> <ycenter> <resolution> <maxiter> <gamma> <filename>");
		System.err.println ("<K> = Number of processors");
		System.err.println ("<width> = Image width (pixels)");
		System.err.println ("<height> = Image height (pixels)");
		System.err.println ("<xcenter> = X coordinate of center point");
		System.err.println ("<ycenter> = Y coordinate of center point");
		System.err.println ("<resolution> = Pixels per unit");
		System.err.println ("<maxiter> = Maximum number of iterations");
		System.err.println ("<gamma> = Used to calculate pixel hues");
		System.err.println ("<filename> = PJG image file name");
		System.exit (1);
		}

	}
