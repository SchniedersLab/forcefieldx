//******************************************************************************
//
// File:    FloydClu.java
// Package: edu.rit.clu.network
// Unit:    Class edu.rit.clu.network.FloydClu
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

import edu.rit.io.DoubleMatrixFile;
import edu.rit.io.Files;

import edu.rit.mp.DoubleBuf;

import edu.rit.pj.Comm;

import edu.rit.util.Range;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/**
 * Class FloydClu is a cluster parallel program that uses Floyd's Algorithm to
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
 * Usage: java -Dpj.np=<I>K</I> edu.rit.clu.network.FloydClu <I>infile</I>
 * <I>outfile</I>
 * <BR><I>K</I> = Number of parallel processes
 * <BR><I>infile</I> = Input distance matrix file
 * <BR><I>outfile</I> = Output distance matrix file
 * <P>
 * The input file (<I>infile</I>) is a binary file written in the format
 * required by class {@linkplain edu.rit.io.DoubleMatrixFile}.
 * <P>
 * Each process writes its own output file containing a row slice of the
 * distance matrix after running Floyd's Algorithm. Each output file is a binary
 * file written in the format required by class {@linkplain
 * edu.rit.io.DoubleMatrixFile}. If <I>outfile</I> is specified as, for example,
 * <TT>"output.dat"</TT>, then the per-process output files are named
 * <TT>"output_0.dat"</TT>, <TT>"output_1.dat"</TT>, and so on through
 * <I>K</I>-1.
 * <P>
 * The computation is performed in parallel in multiple processors. The program
 * measures the total running time (including I/O) and the computation's running
 * time (excluding I/O).
 *
 * @author  Alan Kaminsky
 * @version 07-Jan-2008
 */
public class FloydClu
	{

// Prevent construction.

	private FloydClu()
		{
		}

// Shared variables.

	// World communicator.
	static Comm world;
	static int size;
	static int rank;

	// Number of nodes.
	static int n;

	// Distance matrix.
	static double[][] d;

	// Row broadcast from another process.
	static double[] row_i;
	static DoubleBuf row_i_buf;

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

		// Initialize world communicator.
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();

		// Parse command line arguments.
		if (args.length != 2) usage();
		File infile = new File (args[0]);
		File outfile = new File (args[1]);

		// Prepare to read distance matrix from input file; determine matrix
		// dimensions.
		DoubleMatrixFile in = new DoubleMatrixFile();
		DoubleMatrixFile.Reader reader =
			in.prepareToRead
				(new BufferedInputStream
					(new FileInputStream (infile)));
		d = in.getMatrix();
		n = d.length;

		// Divide distance matrix into equal row slices.
		Range[] ranges = new Range (0, n-1) .subranges (size);
		Range myrange = ranges[rank];
		int mylb = myrange.lb();
		int myub = myrange.ub();

		// Read just this process's row slice of the distance matrix.
		reader.readRowSlice (myrange);
		reader.close();

		// Allocate storage for row broadcast from another process.
		row_i = new double [n];
		row_i_buf = DoubleBuf.buffer (row_i);

		long t2 = System.currentTimeMillis();

		// Run Floyd's Algorithm.
		//     for i = 0 to N-1
		//         for r = 0 to N-1
		//             for c = 0 to N-1
		//                 D[r,c] = min (D[r,c], D[r,i] + D[i,c])
		int i_root = 0;
		for (int i = 0; i < n; ++ i)
			{
			double[] d_i = d[i];

			// Determine which process owns row i.
			if (! ranges[i_root].contains (i)) ++ i_root;

			// Broadcast row i from owner process to all processes.
			if (rank == i_root)
				{
				world.broadcast (i_root, DoubleBuf.buffer (d_i));
				}
			else
				{
				world.broadcast (i_root, row_i_buf);
				d_i = row_i;
				}

			// Inner loops over rows in my slice and over all columns.
			for (int r = mylb; r <= myub; ++ r)
				{
				double[] d_r = d[r];
				for (int c = 0; c < n; ++ c)
					{
					d_r[c] = Math.min (d_r[c], d_r[i] + d_i[c]);
					}
				}
			}

		long t3 = System.currentTimeMillis();

		// Write distance matrix slice to a separate output file in each
		// process.
		DoubleMatrixFile out = new DoubleMatrixFile (n, n, d);
		DoubleMatrixFile.Writer writer =
			out.prepareToWrite
				(new BufferedOutputStream
					(new FileOutputStream
						(Files.fileForRank (outfile, rank))));
		writer.writeRowSlice (myrange);
		writer.close();

		// Stop timing.
		long t4 = System.currentTimeMillis();
		System.out.println ((t2-t1) + " msec pre " + rank);
		System.out.println ((t3-t2) + " msec calc " + rank);
		System.out.println ((t4-t3) + " msec post " + rank);
		System.out.println ((t4-t1) + " msec total " + rank);
		}

// Hidden operations.

	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<K> edu.rit.clu.network.FloydClu <infile> <outfile>");
		System.err.println ("<K> = Number of parallel processes");
		System.err.println ("<infile> = Input distance matrix file");
		System.err.println ("<outfile> = Output distance matrix file");
		System.exit (1);
		}

	}
