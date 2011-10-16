//******************************************************************************
//
// File:    FloydPrint.java
// Package: edu.rit.smp.network
// Unit:    Class edu.rit.smp.network.FloydPrint
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

package edu.rit.smp.network;

import edu.rit.io.DoubleMatrixFile;

import java.io.BufferedInputStream;
import java.io.FileInputStream;

import java.text.DecimalFormat;

/**
 * Class FloydPrint is a main program that prints the contents of a distance
 * matrix input or output file for the {@linkplain FloydSeq}, {@linkplain
 * FloydSmpRow}, and {@linkplain FloydSmpCol} programs.
 * <P>
 * Usage: java edu.rit.smp.network.FloydPrint <I>matrixfile</I>
 * <P>
 * The distance matrix file is a binary file written in the format required by
 * class {@linkplain edu.rit.io.DoubleMatrixFile}.
 *
 * @author  Alan Kaminsky
 * @version 09-Jan-2008
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

		// Read distance matrix file.
		DoubleMatrixFile dmf = new DoubleMatrixFile();
		DoubleMatrixFile.Reader reader =
			dmf.prepareToRead
				(new BufferedInputStream
					(new FileInputStream (matrixfile)));
		reader.read();
		reader.close();

		// Print distance matrix elements.
		int n = dmf.getRowCount();
		double[][] d = dmf.getMatrix();
		for (int r = 0; r < n; ++ r)
			{
			double[] d_r = d[r];
			for (int c = 0; c < n; ++ c)
				{
				StringBuilder buf = new StringBuilder (6);
				buf.append (FMT3.format (d_r[c]));
				while (buf.length() < 6) buf.append (' ');
				System.out.print (buf);
				}
			System.out.println();
			}
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.smp.network.FloydPrint <matrixfile>");
		System.err.println ("<matrixfile> = Distance matrix file");
		System.exit (1);
		}

	}
