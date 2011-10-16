//******************************************************************************
//
// File:    FloydPrint.java
// Package: edu.rit.clu.network
// Unit:    Class edu.rit.clu.network.FloydPrint
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

package edu.rit.clu.network;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;

import java.text.DecimalFormat;

/**
 * Class FloydPrint is a main program that prints the contents of a distance
 * matrix input or output file for the {@linkplain FloydSeq} or {@linkplain
 * FloydClu} programs.
 * <P>
 * Usage: java edu.rit.clu.network.FloydPrint <I>matrixfile</I>
 * <P>
 * The distance matrix file is a binary file written in the format required by
 * class {@linkplain edu.rit.io.DoubleMatrixFile}.
 *
 * @author  Alan Kaminsky
 * @version 07-Jan-2008
 */
public class FloydPrint
	{

// Hidden constants.

	private static final DecimalFormat FMT3 = new DecimalFormat ("0.000");

// Prevent construction.

	private FloydPrint()
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
		if (args.length != 1) usage();
		String matrixfile = args[0];

		// Open distance matrix file.
		DataInputStream in =
			new DataInputStream
				(new BufferedInputStream
					(new FileInputStream (matrixfile)));

		// Read DoubleMatrixFile header information.
		in.readInt(); // Number of matrix rows
		in.readInt(); // Number of matrix columns

		// Read and print each segment until end-of-file.
		try
			{
			for (;;)
				{
				// Read segment header information.
				int RL = in.readInt(); // Segment starting row index
				int CL = in.readInt(); // Segment starting column index
				int M = in.readInt(); // Number of segment rows
				int N = in.readInt(); // Number of segment columns
				System.out.println
					("Rows "+RL+".."+(RL+M-1)+", columns "+CL+".."+(CL+N-1));

				// Read and print distance matrix elements.
				for (int r = 0; r < M; ++ r)
					{
					for (int c = 0; c < N; ++ c)
						{
						StringBuilder buf = new StringBuilder (6);
						buf.append (FMT3.format (in.readDouble()));
						while (buf.length() < 6) buf.append (' ');
						System.out.print (buf);
						}
					System.out.println();
					}
				}
			}
		catch (EOFException exc)
			{
			}

		// Close distance matrix file.
		in.close();
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.clu.network.FloydPrint <matrixfile>");
		System.err.println ("<matrixfile> = Distance matrix file");
		System.exit (1);
		}

	}
