//******************************************************************************
//
// File:    MandelbrotSetGray.java
// Package: edu.rit.smp.fractal
// Unit:    Class edu.rit.smp.fractal.MandelbrotSetGray
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

import edu.rit.image.GrayImageRow;
import edu.rit.image.PJGGrayImage;
import edu.rit.image.PJGImage;

//import edu.rit.pj.Comm;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;

/**
 * Class MandelbrotSetGray is a sequential program that calculates the
 * Mandelbrot Set and creates a grayscale image of it.
 * <P>
 * Usage: java edu.rit.smp.fractal.MandelbrotSetGray <I>width</I> <I>height</I>
 * <I>xcenter</I> <I>ycenter</I> <I>resolution</I> <I>maxiter</I> <I>grayL</I>
 * <I>grayU</I> <I>gamma</I> <I>filename</I>
 * <BR><I>width</I> = Image width (pixels)
 * <BR><I>height</I> = Image height (pixels)
 * <BR><I>xcenter</I> = X coordinate of center point
 * <BR><I>ycenter</I> = Y coordinate of center point
 * <BR><I>resolution</I> = Pixels per unit
 * <BR><I>maxiter</I> = Maximum number of iterations
 * <BR><I>gmin</I> = Used to calculate pixel gray levels
 * <BR><I>gmax</I> = Used to calculate pixel gray levels
 * <BR><I>gslope</I> = Used to calculate pixel gray levels
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
 * <I>c</I> and the pixels' gray values corresponding to the value of <I>i</I>
 * achieved by the iteration. Points in the Mandelbrot set are black, and the
 * other points are shades of gray given by the formula
 * min(<I>gmin</I>+<I>gslope</I>*<I>i</I>/<I>maxiter</I>,<I>gmax</I>), where 0
 * is white and 1 is black. The image is stored in a Parallel Java Graphics
 * (PJG) file specified on the command line.
 * <P>
 * The computation is performed sequentially in a single processor. The program
 * measures the computation's running time, including the time to write the
 * image file. This establishes a benchmark for measuring the computation's
 * running time on a parallel processor.
 *
 * @author  Alan Kaminsky
 * @version 11-Nov-2007
 */
public class MandelbrotSetGray
	{

// Prevent construction.

	private MandelbrotSetGray()
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
	static double gmin;
	static double gmax;
	static double gslope;
	static File filename;

	// Initial pixel offsets from center.
	static int xoffset;
	static int yoffset;

	// Image matrix.
	static byte[][] matrix;
	static PJGGrayImage image;
	static GrayImageRow imagerow;

	// Table of gray levels.
	static float[] graytable;

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
		if (args.length != 10) usage();
		width = Integer.parseInt (args[0]);
		height = Integer.parseInt (args[1]);
		xcenter = Double.parseDouble (args[2]);
		ycenter = Double.parseDouble (args[3]);
		resolution = Double.parseDouble (args[4]);
		maxiter = Integer.parseInt (args[5]);
		gmin = Double.parseDouble (args[6]);
		gmax = Double.parseDouble (args[7]);
		gslope = Double.parseDouble (args[8]);
		filename = new File (args[9]);

		// Initial pixel offsets from center.
		xoffset = -(width - 1) / 2;
		yoffset = (height - 1) / 2;

		// Create image matrix to store results.
		matrix = new byte [height] [width];
		image = new PJGGrayImage (height, width, matrix);
		imagerow = new GrayImageRow (matrix[0]);
		imagerow.setInterpretation (GrayImageRow.ZERO_IS_WHITE);

		// Create table of gray levels for different iteration counts.
		graytable = new float [maxiter+1];
		for (int i = 0; i < maxiter; ++ i)
			{
			graytable[i] = (float) Math.min (gmin + gslope*i/maxiter, gmax);
			}
		graytable[maxiter] = 1.0f;

		long t2 = System.currentTimeMillis();

		// Compute all rows and columns.
		for (int r = 0; r < height; ++ r)
			{
			imagerow.setArray (matrix[r]);
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
				imagerow.setPixel (c, graytable[i]);
				}
			}

		long t3 = System.currentTimeMillis();

		// Write image to file.
		PJGImage.Writer writer =
			image.prepareToWrite
				(new BufferedOutputStream
					(new FileOutputStream (filename)));
		writer.write();
		writer.close();

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
		System.err.println ("Usage: java edu.rit.smp.fractal.MandelbrotSetGray <width> <height> <xcenter> <ycenter> <resolution> <maxiter> <grayL> <grayU> <gamma> <filename>");
		System.err.println ("<width> = Image width (pixels)");
		System.err.println ("<height> = Image height (pixels)");
		System.err.println ("<xcenter> = X coordinate of center point");
		System.err.println ("<ycenter> = Y coordinate of center point");
		System.err.println ("<resolution> = Pixels per unit");
		System.err.println ("<maxiter> = Maximum number of iterations");
		System.err.println ("<grayL> = Used to calculate pixel gray levels");
		System.err.println ("<grayU> = Used to calculate pixel gray levels");
		System.err.println ("<gamma> = Used to calculate pixel gray levels");
		System.err.println ("<filename> = PJG image file name");
		System.exit (1);
		}

	}
