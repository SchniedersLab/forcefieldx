//******************************************************************************
//
// File:    SpinRelaxometryClu.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.SpinRelaxometryClu
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

package edu.rit.mri;

import edu.rit.io.Files;

import edu.rit.mp.ObjectBuf;

import edu.rit.mp.buf.ObjectItemBuf;

import edu.rit.numeric.ArraySeries;
import edu.rit.numeric.Series;

import edu.rit.pj.Comm;
import edu.rit.pj.CommStatus;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;

import edu.rit.util.Range;

import java.io.File;
import java.io.IOException;

import java.util.ArrayList;

/**
 * Class SpinRelaxometryClu is a cluster parallel program that does a spin
 * relaxometry analysis of one or more magnetic resonance images. Each MR
 * image's spin signal data set is stored in a file as defined in class
 * {@linkplain SignalDataSetWriter}. The program does the following for each
 * spin signal data set file. Using an instance of class {@linkplain
 * SignalDataSetReader}, the program reads the spin signal data for each pixel
 * in the image. Using class {@linkplain PixelAnalysis}, the program does a spin
 * relaxometry analysis on each pixel and computes the tissues for each pixel.
 * Using an instance of class {@linkplain TissuesDataSetWriter}, the program
 * writes the tissues data into another file. For example, if the input spin
 * signal data set files are named <TT>image1.dat</TT>, <TT>image2.dat</TT>, and
 * so on, the output tissues data set files are named
 * <TT>tissues_image1.dat</TT>, <TT>tissues_image2.dat</TT>, and so on.
 * <P>
 * The program uses the master-worker pattern for load balancing. The program
 * uses the parallel input files pattern for reduced message passing. The master
 * partitions each image into chunks of 100 pixels and sends the chunks to the
 * workers. In parallel, the workers read the input spin signal data set file,
 * do the spin relaxometry analysis calculations for each pixel in the chunk,
 * and send the calculated tissues data back to the master. The master writes
 * the tissues data into the output tissues data set file.
 * <P>
 * Usage: java -Dpj.np=<I>K</I> edu.rit.mri.SpinRelaxometryClu <I>R1_lower</I>
 * <I>R1_upper</I> <I>N</I> <I>signalfile</I> [ <I>signalfile</I> . . . ]
 * <BR><I>K</I> = Number of parallel processes
 * <BR><I>R1_lower</I> = Lower <I>R</I>1 spin-lattice relaxation rate (1/sec)
 * <BR><I>R1_upper</I> = Upper <I>R</I>1 spin-lattice relaxation rate (1/sec)
 * <BR><I>N</I> = Number of <I>R</I>1 intervals
 * <BR><I>signalfile</I> = Input spin signal data set file
 *
 * @author  Alan Kaminsky
 * @version 25-Jun-2008
 */
public class SpinRelaxometryClu
	{

// Global variables.

	// World communicator.
	static Comm world;
	static int size;
	static int rank;

	// Command line arguments.
	static double R1_lower;
	static double R1_upper;
	static int N;
	static String[] signalfilename;

// Prevent construction.

	private SpinRelaxometryClu()
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
		// Start timing.
		long t1 = System.currentTimeMillis();

		// Initialize world communicator.
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();

		// Parse command line arguments.
		if (args.length < 4) usage();
		R1_lower = Double.parseDouble (args[0]);
		R1_upper = Double.parseDouble (args[1]);
		N = Integer.parseInt (args[2]);
		signalfilename = new String [args.length-3];
		System.arraycopy (args, 3, signalfilename, 0, args.length-3);

		// Process 0 executes the master section and the worker section in
		// parallel threads.
		if (rank == 0)
			{
			new ParallelTeam(2).execute (new ParallelRegion()
				{
				public void run() throws Exception
					{
					execute (new ParallelSection()
						{
						public void run() throws Exception
							{
							masterSection();
							}
						},
					new ParallelSection()
						{
						public void run() throws Exception
							{
							workerSection();
							}
						});
					}
				});
			}

		// Processes 1 and up execute just the worker section.
		else
			{
			workerSection();
			}

		// Stop timing.
		long t2 = System.currentTimeMillis();
		System.out.println ((t2-t1)+" msec "+rank);
		}

// Hidden operations.

	/**
	 * Execute the master section.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private static void masterSection()
		throws IOException
		{
		int worker;
		PixelChunk chunk;

		// Set up array of writers for output tissues data sets.
		TissuesDataSetWriter[] writer =
			new TissuesDataSetWriter [signalfilename.length];
		for (int f = 0; f < writer.length; ++ f)
			{
			SignalDataSetReader reader =
				new SignalDataSetReader (new File (signalfilename[f]));
			int H = reader.getHeight();
			int W = reader.getWidth();
			reader.close();
			writer[f] =
				new TissuesDataSetWriter
					(new File
						(Files.fileNamePrepend (signalfilename[f], "tissues_")),
					 H, W);
			}

		// Set up schedule to analyze pixels in chunks of 100.
		PixelSchedule schedule = new PixelSchedule (100, signalfilename);

		// Send initial chunk to each worker. If null, no more work for that
		// worker. Keep count of active workers.
		int activeWorkers = size;
		for (worker = 0; worker < size; ++ worker)
			{
			chunk = schedule.next();
			world.send (worker, ObjectBuf.buffer (chunk));
			if (chunk == null) -- activeWorkers;
			}

		// Repeat until all workers have finished.
		while (activeWorkers > 0)
			{
			// Receive a chunk of pixel tissues from any worker.
			ObjectItemBuf<PixelTissues[]> buf = ObjectBuf.buffer();
			CommStatus status = world.receive (null, buf);
			worker = status.fromRank;

			// Send next chunk to that specific worker. If null, no more work.
			chunk = schedule.next();
			world.send (worker, ObjectBuf.buffer (chunk));
			if (chunk == null) -- activeWorkers;

			// Record pixel tissues in output tissues data set.
			for (PixelTissues tissues : buf.item)
				{
				if (tissues != null)
					{
					writer[tissues.fileIndex()].addPixelTissues (tissues);
					}
				}
			}

		// All done.
		for (int f = 0; f < writer.length; ++ f)
			{
			writer[f].close();
			}
		}

	/**
	 * Execute the worker section.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private static void workerSection()
		throws IOException
		{
		int fileIndex = -1;
		SignalDataSetReader reader = null;
		Series t_series = null;
		int M = 0;
		double[][] A = null;

		// Compute spin relaxation rates.
		double[] R1 = new double [N+1];
		double log_R1_lower = Math.log (R1_lower);
		double log_R1_upper = Math.log (R1_upper);
		double interval = (log_R1_upper - log_R1_lower)/N;
		for (int j = 0; j <= N; ++ j)
			{
			R1[j] = Math.exp (log_R1_lower + j*interval);
			}
		ArraySeries R1_series = new ArraySeries (R1);

		// Set up lists to receive analysis results.
		ArrayList<Double> rho_list = new ArrayList<Double>();
		ArrayList<Double> R1_list = new ArrayList<Double>();

		// Repeat until no more work.
		workerloop: for (;;)
			{
			// Receive a chunk of pixel indexes from the master. If null, no
			// more work.
			ObjectItemBuf<PixelChunk> buf = ObjectBuf.buffer();
			world.receive (0, buf);
			PixelChunk chunk = buf.item;
			if (chunk == null) break workerloop;
			int f = chunk.fileIndex();
			int lb = chunk.pixelIndex();
			int len = chunk.pixelCount();

			// If we are now working on a different file:
			if (f != fileIndex)
				{
				// Close old file.
				if (reader != null) reader.close();

				// Open new file.
				fileIndex = f;
				reader = new SignalDataSetReader (new File (signalfilename[f]));

				// Get time series.
				t_series = reader.getTimeSeries();
				M = t_series.length();

				// Compute design matrix.
				A = new double [M] [N+1];
				for (int i = 0; i < M; ++ i)
					{
					double[] A_i = A[i];
					double t_i = t_series.x(i);
					for (int j = 0; j <= N; ++ j)
						{
						A_i[j] = SpinSignal.S (R1[j], t_i);
						}
					}
				}

			// Set up array of pixel tissues to hold analysis results.
			PixelTissues[] tissues = new PixelTissues [len];

			// Process all pixels in chunk.
			for (int i = 0; i < len; ++ i)
				{
				int index = lb + i;
				PixelSignal signal_i = reader.getPixelSignal (index);
				if (signal_i != null)
					{
					// Get measured spin signal.
					Series S_series = signal_i.S_measured();

					// Do the spin relaxometry analysis.
					PixelAnalysis.analyze
						(t_series, S_series, R1_series, A, rho_list, R1_list);

					// Record analysis results.
					tissues[i] = new PixelTissues (f, index, rho_list, R1_list);
					}
				}

			// Send chunk of pixel tissues to the master.
			world.send (0, ObjectBuf.objectBuffer (tissues));
			}

		// All done.
		if (reader != null) reader.close();
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<K> edu.rit.mri.SpinRelaxometryClu <R1_lower> <R1_upper> <N> <signalfile> [<signalfile> ...]");
		System.err.println ("<K> = Number of parallel processes");
		System.err.println ("<R1_lower> = Lower R1 spin-lattice relaxation rate (1/sec)");
		System.err.println ("<R1_upper> = Upper R1 spin-lattice relaxation rate (1/sec)");
		System.err.println ("<N> = Number of R1 intervals");
		System.err.println ("<signalfile> = Input spin signal data set file");
		System.exit (1);
		}

	}
