//******************************************************************************
//
// File:    FloydSeq.java
// Package: edu.rit.hyb.network
// Unit:    Class edu.rit.hyb.network.FloydSeq
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

package edu.rit.hyb.network;

import edu.rit.io.DoubleMatrixFile;

import edu.rit.pj.Comm;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/**
 * Class FloydSeq is a sequential program that uses Floyd's Algorithm to
 * calculate the length of the shortest path from each node to every other node
 * in a network, given the distance from each node to its adjacent nodes.
 * <P>
 * Floyd's Algorithm's running time is <I>O</I>(<I>N</I><SUP>3</SUP>), where
 * <I>N</I> is the number of nodes. The algorithm is as follows. On input,
 * <I>D</I> is an <I>N</I>x<I>N</I> matrix where <I>D[i,j]</I> is the distance
 * from node <I>i</I> to adjacent node <I>j</I>; if node <I>j</I> is not
 * adjacent to node <I>i</I>, then <I>D[i,j]</I> is infinity. On output,
 * <I>D[i,j]</I> has been replaced by the length of the shortest path from node
 * <I>i</I> to node <I>j</I>; if there is no path from node <I>i</I> to node
 * <I>j</I>, then <I>D[i,j]</I> is infinity.
 * <PRE>
 *     for i = 0 to N-1
 *         for r = 0 to N-1
 *             for c = 0 to N-1
 *                 D[r,c] = min (D[r,c], D[r,i] + D[i,c])
 * </PRE>
 * <P>
 * Usage: java edu.rit.hyb.network.FloydSeq <I>infile</I> <I>outfile</I>
 * <BR><I>infile</I> = Input distance matrix file
 * <BR><I>outfile</I> = Output distance matrix file
 * <P>
 * The input file (<I>infile</I>) is a binary file written in the format
 * required by class {@linkplain edu.rit.io.DoubleMatrixFile}.
 * <P>
 * The output file (<I>outfile</I>) is a binary file written in the format
 * required by class {@linkplain edu.rit.io.DoubleMatrixFile} containing the
 * distance matrix after running Floyd's Algorithm.
 * <P>
 * The computation is performed sequentially in a single processor. The program
 * measures the total running time (including I/O) and the computation's running
 * time (excluding I/O). This establishes a benchmark for measuring the running
 * time on a parallel processor.
 *
 * @author  Alan Kaminsky
 * @version 30-May-2008
 */
public class FloydSeq
	{

// Prevent construction.

	private FloydSeq()
		{
		}

// Shared variables.

	// Number of nodes.
	static int n;

	// Distance matrix.
	static double[][] d;

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Throwable
		{
		// Start timing.
		long t1 = System.currentTimeMillis();

		// Initialize middleware.
		Comm.init (args);

		// Parse command line arguments.
		if (args.length != 2) usage();
		File infile = new File (args[0]);
		File outfile = new File (args[1]);

		// Read distance matrix from input file.
		DoubleMatrixFile in = new DoubleMatrixFile();
		DoubleMatrixFile.Reader reader =
			in.prepareToRead
				(new BufferedInputStream
					(new FileInputStream (infile)));
		reader.read();
		reader.close();
		d = in.getMatrix();
		n = d.length;

		// Run Floyd's Algorithm.
		//     for i = 0 to N-1
		//         for r = 0 to N-1
		//             for c = 0 to N-1
		//                 D[r,c] = min (D[r,c], D[r,i] + D[i,c])
		long t2 = System.currentTimeMillis();
		for (int i = 0; i < n; ++ i)
			{
			double[] d_i = d[i];
			for (int r = 0; r < n; ++ r)
				{
				double[] d_r = d[r];
				for (int c = 0; c < n; ++ c)
					{
					d_r[c] = Math.min (d_r[c], d_r[i] + d_i[c]);
					}
				}
			}
		long t3 = System.currentTimeMillis();

		// Write distance matrix to output file.
		DoubleMatrixFile out = new DoubleMatrixFile (n, n, d);
		DoubleMatrixFile.Writer writer =
			out.prepareToWrite
				(new BufferedOutputStream
					(new FileOutputStream (outfile)));
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

	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.hyb.network.FloydSeq <infile> <outfile>");
		System.err.println ("<infile> = Input distance matrix file");
		System.err.println ("<outfile> = Output distance matrix file");
		System.exit (1);
		}

	}
