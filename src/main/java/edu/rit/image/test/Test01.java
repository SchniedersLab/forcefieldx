//******************************************************************************
//
// File:    Test01.java
// Package: edu.rit.image.test
// Unit:    Class edu.rit.image.test.Test01
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

package edu.rit.image.test;

import edu.rit.color.HSB;

import edu.rit.image.PJGColorImage;
import edu.rit.image.PJGImage;

import edu.rit.util.Hex;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/**
 * Class Test01 is a unit test program for class {@linkplain
 * edu.rit.image.PJGColorImage}. The program calculates an image of the
 * Mandelbrot Set and stores the image in a Parallel Java Graphics (PJG) file.
 * Then the program reads the PJG file back in and tests whether the pixel data
 * that came in is the same as the pixel data that went out. The image is
 * <I>width</I> pixels wide and <I>width</I> pixels high.
 * <P>
 * Usage: java edu.rit.image.test.Test01 <I>width</I> <I>filename</I>
 * <BR><I>width</I> = Image width and height (pixels)
 * <BR><I>filename</I> = PJG image file name
 *
 * @author  Alan Kaminsky
 * @version 02-Nov-2007
 */
public class Test01
	{

// Prevent construction.

	private Test01()
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

	// Image matrices.
	static int[][] matrix;
	static int[][] matrix2;

	// Table of hues.
	static int[] huetable;

// Main program.

	/**
	 * Unit test main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		long t1 = System.currentTimeMillis();

		// Validate command line arguments.
		if (args.length != 2) usage();
		width = Integer.parseInt (args[0]);
		height = width;
		xcenter = -0.75;
		ycenter = 0.0;
		resolution = width * 0.375;
		maxiter = 1000;
		gamma = 0.4;
		filename = new File (args[1]);

		// Initial pixel offsets from center.
		xoffset = -(width - 1) / 2;
		yoffset = (height - 1) / 2;

		// Create image matrix to store results.
		matrix = new int [height] [width];

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
		System.out.println ((t2-t1) + " msec -- Calculating image");

		// Compute all rows and columns.
		for (int r = 0; r < height; ++ r)
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
		System.out.println ((t3-t1) + " msec -- Writing PJG file");

		// Write image to file.
		PJGColorImage image = new PJGColorImage (height, width, matrix);
		PJGImage.Writer writer =
			image.prepareToWrite
				(new BufferedOutputStream
					(new FileOutputStream (filename)));
		writer.write();
		writer.close();

		long t4 = System.currentTimeMillis();
		System.out.println ((t4-t1) + " msec -- Reading PJG file");

		// Read image from file into a different matrix.
		PJGColorImage image2 = new PJGColorImage();
		PJGImage.Reader reader =
			image2.prepareToRead
				(new BufferedInputStream
					(new FileInputStream (filename)));
		reader.read();
		reader.close();

		long t5 = System.currentTimeMillis();
		System.out.println ((t5-t1) + " msec -- Comparing images");

		// Compare output matrix with input matrix.
		matrix2 = image2.getMatrix();
		for (int r = 0; r < width; ++ r)
			{
			int[] matrix_r = matrix[r];
			int[] matrix2_r = matrix2[r];
			for (int c = 0; c < width; ++ c)
				{
				if (matrix_r[c] != matrix2_r[c])
					{
					System.out.print ("matrix[");
					System.out.print (r);
					System.out.print ("][");
					System.out.print (c);
					System.out.print ("] = ");
					System.out.print (Hex.toString (matrix_r[c]));
					System.out.print (", matrix2[");
					System.out.print (r);
					System.out.print ("][");
					System.out.print (c);
					System.out.print ("] = ");
					System.out.print (Hex.toString (matrix2_r[c]));
					System.out.println();
					}
				}
			}

		long t6 = System.currentTimeMillis();
		System.out.println ((t6-t1) + " msec -- Done");
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.image.test.Test01 <width> <filename>");
		System.err.println ("<width> = Image width and height (pixels)");
		System.err.println ("<filename> = PNG image file name");
		System.exit (1);
		}

	}
