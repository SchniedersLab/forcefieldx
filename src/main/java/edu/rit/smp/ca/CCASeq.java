//******************************************************************************
//
// File:    CCASeq.java
// Package: edu.rit.smp.ca
// Unit:    Class edu.rit.smp.ca.CCASeq
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

package edu.rit.smp.ca;

import edu.rit.image.GrayImageRow;
import edu.rit.image.PJGGrayImage;
import edu.rit.image.PJGImage;

import edu.rit.numeric.BigRational;

import edu.rit.util.Range;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Class CCASeq is a sequential program that calculates the evolution of a
 * continuous cellular automaton and stores the result in a grayscale PJG image
 * file.
 * <P>
 * Usage: java edu.rit.smp.ca.CCASeq <I>C</I> <I>S</I> <I>A</I> <I>B</I>
 * <I>imagefile</I>
 * <BR><I>C</I> = Number of cells (&gt;= 1)
 * <BR><I>S</I> = Number of time steps (&gt;= 1)
 * <BR><I>A</I> = Multiplicand in update formula (rational number)
 * <BR><I>B</I> = Addend in update formula (rational number)
 * <BR><I>imagefile</I> = Output PJG image file name
 * <P>
 * The cellular automaton (CA) consists of an array of <I>C</I> cells. Each
 * cell's value is a rational number in the range 0 to 1. A cell's next value is
 * computed as follows: Compute the average of the cell itself, the cell's left
 * neighbor, and the cell's right neighbor (wraparound boundary conditions);
 * multiply by <I>A</I>; add <I>B</I>; and take the fractional part, yielding a
 * result in the range 0 to 1. The CA's initial state is all cells 0, except the
 * middle cell is 1. The program evolves the CA's initial state for <I>S</I>
 * time steps and generates a grayscale image with <I>S</I>+1 rows and <I>C</I>
 * columns. The first row of the image corresponds to the CA's initial state,
 * and each subsequent row corresponds to the CA's state after the next time
 * step. In the image, each pixel's gray value is proportional to the
 * corresponding cell's value, with 0 being white and 1 being black. The image
 * is stored in a Parallel Java Graphics (PJG) file. For further information
 * about the PJG format, see class {@linkplain edu.rit.image.PJGImage}.
 *
 *
 * @author  Alan Kaminsky
 * @version 06-Jul-2009
 */
public class CCASeq
	{

// Prevent construction.

	private CCASeq()
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
	static File imagefile;

	// Old and new cell arrays.
	static BigRational[] currentCell;
	static BigRational[] nextCell;

	// Grayscale image matrix.
	static byte[][] pixelmatrix;
	static PJGGrayImage image;
	static PJGImage.Writer writer;

	// One row of the grayscale image matrix.
	static byte[] pixelrow;
	static GrayImageRow imagerow;

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Start timing.
		long t1 = System.currentTimeMillis();

		// Parse command line arguments.
		if (args.length != 5) usage();
		C = Integer.parseInt (args[0]);
		S = Integer.parseInt (args[1]);
		A = new BigRational (args[2]) .mul (ONE_THIRD);
		B = new BigRational (args[3]);
		imagefile = new File (args[4]);

		// Allocate storage for old and new cell arrays. Initialize all cells to
		// 0, except center cell to 1.
		currentCell = new BigRational [C];
		nextCell = new BigRational [C];
		for (int i = 0; i < C; ++ i)
			{
			currentCell[i] = new BigRational();
			nextCell[i] = new BigRational();
			}
		currentCell[C/2].assign (ONE);

		// Set up pixel matrix, image, and image writer.
		pixelmatrix = new byte [S+1] [];
		image = new PJGGrayImage (S+1, C, pixelmatrix);
		writer =
			image.prepareToWrite
				(new BufferedOutputStream
					(new FileOutputStream (imagefile)));

		// Allocate storage for one pixel matrix row.
		pixelrow = new byte [C];
		imagerow = new GrayImageRow (pixelrow);
		imagerow.setInterpretation (PJGGrayImage.ZERO_IS_WHITE);

		// Do S time steps.
		for (int s = 0; s < S; ++ s)
			{
			// Calculate next state of each cell.
			for (int i = 0; i < C; ++ i)
				{
				nextCell[i]
					.assign (currentCell[i])
					.add (currentCell[(i-1+C)%C])
					.add (currentCell[(i+1)%C])
					.mul (A)
					.add (B)
					.normalize()
					.fracPart();
				}

			// Write current CA state to image file.
			writeCurrentCell (s);

			// Advance one time step -- swap old and new cell arrays.
			BigRational[] tmp = currentCell;
			currentCell = nextCell;
			nextCell = tmp;
			}

		// Write final CA state to image file.
		writeCurrentCell (S);
		writer.close();

		// Stop timing.
		long t2 = System.currentTimeMillis();
		System.out.println ((t2-t1) + " msec total");
		}

// Hidden operations.

	/**
	 * Write the current cell array to the given row of the image file.
	 *
	 * @param  r  Row index.
	 */
	private static void writeCurrentCell
		(int r)
		throws IOException
		{
		// Set image row's gray values based on current cell states.
		for (int i = 0; i < C; ++ i)
			{
			imagerow.setPixel (i, currentCell[i].floatValue());
			}

		// Set row r of the pixel matrix.
		pixelmatrix[r] = pixelrow;

		// Write row-r slice of the image to the image file.
		writer.writeRowSlice (new Range (r, r));
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.smp.ca.CCASeq <C> <S> <A> <B> <imagefile>");
		System.err.println ("<C> = Number of cells (>= 1)");
		System.err.println ("<S> = Number of time steps (>= 1)");
		System.err.println ("<A> = Multiplicand in update formula (rational number)");
		System.err.println ("<B> = Addend in update formula (rational number)");
		System.err.println ("<imagefile> = Output PJG image file name");
		System.exit (1);
		}

	}
