//******************************************************************************
//
// File:    Test03.java
// Package: edu.rit.image.test
// Unit:    Class edu.rit.image.test.Test03
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

import edu.rit.image.PJGGrayImage;
import edu.rit.image.PJGImage;

import edu.rit.util.Hex;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import java.util.Arrays;

/**
 * Class Test03 is a unit test program for class {@linkplain
 * edu.rit.image.PJGGrayImage}. The program calculates an image of a 2-D
 * Gaussian function and stores the image in a Parallel Java Graphics (PJG)
 * file. Then the program reads the PJG file back in and tests whether the pixel
 * data that came in is the same as the pixel data that went out. The image is
 * <I>width</I> pixels high and <I>width</I> pixels wide.
 * <P>
 * Usage: java edu.rit.image.test.Test03 <I>width</I> <I>filename</I>
 * <BR><I>width</I> = Image height and width (pixels)
 * <BR><I>filename</I> = PJG image file name
 *
 * @author  Alan Kaminsky
 * @version 10-Nov-2007
 */
public class Test03
	{

// Prevent construction.

	private Test03()
		{
		}

// Shared variables.

	// Command line arguments.
	static int width;
	static File filename;

	// Output image matrix.
	static byte[][] matrix;
	static PJGGrayImage image;

	// Input image matrix.
	static byte[][] matrix2;
	static PJGGrayImage image2;

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Parse command line arguments.
		if (args.length != 2) usage();
		width = Integer.parseInt (args[0]);
		filename = new File (args[1]);

		// Allocate storage for output image.
		matrix = new byte [width] [width];
		image = new PJGGrayImage (width, width, matrix);
		image.setInterpretation (PJGGrayImage.ZERO_IS_BLACK);

		// Calculate image of 2-D Gaussian.
		for (int r = 0; r < width; ++ r)
			{
			double y = ((double) r) / ((double) width) * 4.0 - 2.0;
			for (int c = 0; c < width; ++ c)
				{
				double x = ((double) c) / ((double) width) * 4.0 - 2.0;
				double z = Math.exp (- x*x - y*y);
				image.setPixel (r, c, (float) z);
				}
			}

		// Write output image to PJG file.
		PJGImage.Writer writer =
			image.prepareToWrite
				(new BufferedOutputStream
					(new FileOutputStream (filename)));
		writer.write();
		writer.close();

		// Read image from file into a different matrix.
		image2 = new PJGGrayImage();
		PJGImage.Reader reader =
			image2.prepareToRead
				(new BufferedInputStream
					(new FileInputStream (filename)));
		reader.read();
		reader.close();

		// Compare output matrix with input matrix.
		matrix2 = image2.getMatrix();
		for (int r = 0; r < width; ++ r)
			{
			byte[] matrix_r = matrix[r];
			byte[] matrix2_r = matrix2[r];
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
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.image.test.Test03 <width> <filename>");
		System.err.println ("<width> = Image height and width (pixels)");
		System.err.println ("<filename> = PJG image file name");
		System.exit (1);
		}

	}
