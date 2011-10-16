//******************************************************************************
//
// File:    Test02.java
// Package: edu.rit.image.test
// Unit:    Class edu.rit.image.test.Test02
//
// This Java source file is copyright (C) 2009 by Alan Kaminsky. All rights
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

import edu.rit.numeric.BigRational;

import edu.rit.util.Hex;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Class Test02 is a unit test program for class {@linkplain
 * edu.rit.image.PJGGrayImage}. The program calculates an image of the evolution
 * of a continuous cellular automaton and stores the image in a Parallel Java
 * Graphics (PJG) file. Then the program reads the PJG file back in and tests
 * whether the pixel data that came in is the same as the pixel data that went
 * out. The image is <I>height</I> pixels high and <I>width</I> pixels wide.
 * <P>
 * Usage: java edu.rit.image.test.Test02 <I>height</I> <I>width</I>
 * <I>filename</I>
 * <BR><I>height</I> = Image height (pixels)
 * <BR><I>width</I> = Image width (pixels)
 * <BR><I>filename</I> = PJG image file name
 *
 * @author  Alan Kaminsky
 * @version 06-Jul-2009
 */
public class Test02
	{

// Prevent construction.

	private Test02()
		{
		}

// Shared variables.

	// Constants.
	static final BigRational ZERO = new BigRational ("0");
	static final BigRational ONE = new BigRational ("1");
	static final BigRational ONE_THIRD = new BigRational ("1/3");

	// Command line arguments.
	static int C;
	static int S;
	static BigRational A;
	static BigRational B;
	static File filename;

	// Old and new cell arrays.
	static BigRational[] currCell;
	static BigRational[] nextCell;

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
		if (args.length != 3) usage();
		C = Integer.parseInt (args[1]);
		S = Integer.parseInt (args[0]) - 1;
		A = ONE_THIRD;
		B = new BigRational ("9/10");
		filename = new File (args[2]);

		// Allocate storage for old and new cell arrays. Initialize all cells to
		// 0, except center cell to 1.
		currCell = new BigRational [C];
		nextCell = new BigRational [C];
		for (int i = 0; i < C; ++ i)
			{
			currCell[i] = new BigRational();
			nextCell[i] = new BigRational();
			}
		currCell[C/2].assign (ONE);

		// Allocate storage for output image.
		matrix = new byte [S+1] [C];
		image = new PJGGrayImage (S+1, C, matrix);
		image.setInterpretation (PJGGrayImage.ZERO_IS_WHITE);

		// Do S time steps.
		for (int s = 0; s < S; ++ s)
			{
			// Calculate next state of each cell.
			for (int i = 0; i < C; ++ i)
				{
				nextCell[i]
					.assign (currCell[i])
					.add (currCell[(i-1+C)%C])
					.add (currCell[(i+1)%C])
					.mul (A)
					.add (B)
					.normalize()
					.fracPart();
				}

			// Write current CA state to image matrix.
			writeCurrCell (s);

			// Advance one time step -- swap old and new cell arrays.
			BigRational[] tmp = currCell;
			currCell = nextCell;
			nextCell = tmp;
			}

		// Write final CA state to image matrix.
		writeCurrCell (S);

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
		for (int r = 0; r <= S; ++ r)
			{
			byte[] matrix_r = matrix[r];
			byte[] matrix2_r = matrix2[r];
			for (int c = 0; c < C; ++ c)
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
	 * Write the current cell array to the image file.
	 *
	 * @param  s  Step (image row) index.
	 */
	private static void writeCurrCell
		(int s)
		throws IOException
		{
		for (int i = 0; i < C; ++ i)
			{
			image.setPixel (s, i, currCell[i].floatValue());
			}
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.image.test.Test02 <height> <width> <filename>");
		System.err.println ("<height> = Image height (pixels)");
		System.err.println ("<width> = Image width (pixels)");
		System.err.println ("<filename> = PJG image file name");
		System.exit (1);
		}

	}
