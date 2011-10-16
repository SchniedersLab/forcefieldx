//******************************************************************************
//
// File:    MandelbrotSetClu.java
// Package: edu.rit.clu.fractal
// Unit:    Class edu.rit.clu.fractal.MandelbrotSetClu
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

package edu.rit.clu.fractal;

import edu.rit.color.HSB;

import edu.rit.image.PJGColorImage;
import edu.rit.image.PJGImage;

import edu.rit.mp.IntegerBuf;

import edu.rit.pj.Comm;

import edu.rit.util.Arrays;
import edu.rit.util.Range;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;

/**
 * Class MandelbrotSetClu is a cluster parallel program that calculates the
 * Mandelbrot Set. Each process in the program calculates a fixed row slice of
 * the Mandelbrot Set image. The slices are gathered into process 0, which then
 * writes the image to a file.
 * <P>
 * Usage: java -Dpj.np=<I>K</I> edu.rit.clu.fractal.MandelbrotSetClu
 * <I>width</I> <I>height</I> <I>xcenter</I> <I>ycenter</I> <I>resolution</I>
 * <I>maxiter</I> <I>gamma</I> <I>filename</I>
 * <BR><I>K</I> = Number of parallel processes
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
 * image file.
 *
 * @author  Alan Kaminsky
 * @version 02-Nov-2007
 */
public class MandelbrotSetClu
	{

// Prevent construction.

	private MandelbrotSetClu()
		{
		}

// Program shared variables.

	// Communicator.
	static Comm world;
	static int size;
	static int rank;

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
	static Range[] ranges;
	static Range myrange;
	static int mylb;
	static int myub;

	// Communication buffers.
	static IntegerBuf[] slices;
	static IntegerBuf myslice;

	// Table of hues.
	static int[] huetable;

// Main program.

	/**
	 * Mandelbrot Set main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Start timing.
		long t1 = System.currentTimeMillis();

		// Initialize middleware.
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();

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

		// Create image matrix to store results; the full matrix in process 0,
		// one row slice of the matrix in the other processes.
		matrix = new int [height] [];
		ranges = new Range (0, height-1) .subranges (size);
		myrange = ranges[rank];
		mylb = myrange.lb();
		myub = myrange.ub();
		if (rank == 0)
			{
			Arrays.allocate (matrix, width);
			}
		else
			{
			Arrays.allocate (matrix, myrange, width);
			}

		// Set up communication buffers.
		slices = IntegerBuf.rowSliceBuffers (matrix, ranges);
		myslice = slices[rank];

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

		long t2 = System.currentTimeMillis();

		// Compute all rows and columns.
		for (int r = mylb; r <= myub; ++ r)
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

		long t3 = System.currentTimeMillis();

		// Gather all matrix row slices into process 0.
		world.gather (0, myslice, slices);

		// Write image to PJG file in process 0.
		if (rank == 0)
			{
			image = new PJGColorImage (height, width, matrix);
			PJGImage.Writer writer =
				image.prepareToWrite
					(new BufferedOutputStream
						(new FileOutputStream (filename)));
			writer.write();
			writer.close();
			}

		// Stop timing.
		long t4 = System.currentTimeMillis();
		System.out.println ((t2-t1) + " msec pre " + rank);
		System.out.println ((t3-t2) + " msec calc " + rank);
		System.out.println ((t4-t3) + " msec post " + rank);
		System.out.println ((t4-t1) + " msec total " + rank);
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<K> edu.rit.clu.fractal.MandelbrotSetClu <width> <height> <xcenter> <ycenter> <resolution> <maxiter> <gamma> <filename>");
		System.err.println ("<K> = Number of parallel processes");
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
