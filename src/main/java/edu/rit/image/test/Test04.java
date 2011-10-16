//******************************************************************************
//
// File:    Test04.java
// Package: edu.rit.image.test
// Unit:    Class edu.rit.image.test.Test04
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

package edu.rit.image.test;

import edu.rit.image.PJGHueImage;
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
 * Class Test04 is a unit test program for class {@linkplain
 * edu.rit.image.PJGHueImage}. The program calculates an image with hues varying
 * continuously along the horizontal axis and brightnesses varying continuously
 * along the vertical axis and stores the image in a Parallel Java Graphics
 * (PJG) file. Then the program reads the PJG file back in and tests whether the
 * pixel data that came in is the same as the pixel data that went out. The
 * image is <I>width</I> pixels high and <I>width</I> pixels wide.
 * <P>
 * Usage: java edu.rit.image.test.Test04 <I>width</I> <I>filename</I>
 * <BR><I>width</I> = Image height and width (pixels)
 * <BR><I>filename</I> = PJG image file name
 *
 * @author  Alan Kaminsky
 * @version 08-Apr-2008
 */
public class Test04
	{

// Prevent construction.

	private Test04()
		{
		}

// Shared variables.

	// Command line arguments.
	static int width;
	static File filename;

	// Output image matrix.
	static int[][] matrix;
	static PJGHueImage image;

	// Input image matrix.
	static int[][] matrix2;
	static PJGHueImage image2;

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
		matrix = new int [width] [width];
		image = new PJGHueImage (width, width, matrix);

		// Calculate image with continuously varying brightnesses and hues.
		for (int r = 0; r < width; ++ r)
			{
			float bri = 1.0f - ((float) r) / ((float) width);
			for (int c = 0; c < width; ++ c)
				{
				float hue = ((float) c) / ((float) width);
				image.setPixelHSB (r, c, hue, 1.0f, bri);
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
		image2 = new PJGHueImage();
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
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.image.test.Test04 <width> <filename>");
		System.err.println ("<width> = Image height and width (pixels)");
		System.err.println ("<filename> = PJG image file name");
		System.exit (1);
		}

	}
