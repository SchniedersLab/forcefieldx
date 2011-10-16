//******************************************************************************
//
// File:    SpinRelaxometrySeq.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.SpinRelaxometrySeq
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

import edu.rit.numeric.ArraySeries;
import edu.rit.numeric.Series;

import edu.rit.pj.Comm;

import java.io.File;

import java.util.ArrayList;

/**
 * Class SpinRelaxometrySeq is a sequential program that does a spin relaxometry
 * analysis of one or more magnetic resonance images. Each MR image's spin
 * signal data set is stored in a file as defined in class {@linkplain
 * SignalDataSetWriter}. The program does the following for each spin signal
 * data set file. Using an instance of class {@linkplain SignalDataSetReader},
 * the program reads the spin signal data for each pixel in the image. Using
 * class {@linkplain PixelAnalysis}, the program does a spin relaxometry
 * analysis on each pixel and computes the tissues for each pixel. Using an
 * instance of class {@linkplain TissuesDataSetWriter}, the program writes the
 * tissues data into another file. For example, if the input spin signal data
 * set files are named <TT>image1.dat</TT>, <TT>image2.dat</TT>, and so on, the
 * output tissues data set files are named <TT>tissues_image1.dat</TT>,
 * <TT>tissues_image2.dat</TT>, and so on.
 * <P>
 * Usage: java edu.rit.mri.SpinRelaxometrySeq <I>R1_lower</I> <I>R1_upper</I>
 * <I>N</I> <I>signalfile</I> [ <I>signalfile</I> . . . ]
 * <BR><I>R1_lower</I> = Lower <I>R</I>1 spin-lattice relaxation rate (1/sec)
 * <BR><I>R1_upper</I> = Upper <I>R</I>1 spin-lattice relaxation rate (1/sec)
 * <BR><I>N</I> = Number of <I>R</I>1 intervals
 * <BR><I>signalfile</I> = Input spin signal data set file
 *
 * @author  Alan Kaminsky
 * @version 25-Jun-2008
 */
public class SpinRelaxometrySeq
	{

// Prevent construction.

	private SpinRelaxometrySeq()
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

		Comm.init (args);

		// Parse command line arguments.
		if (args.length < 4) usage();
		double R1_lower = Double.parseDouble (args[0]);
		double R1_upper = Double.parseDouble (args[1]);
		int N = Integer.parseInt (args[2]);
		String[] signalfilename = new String [args.length-3];
		System.arraycopy (args, 3, signalfilename, 0, args.length-3);

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

		// Analyze each input spin signal data set file.
		for (int f = 0; f < signalfilename.length; ++ f)
			{
			File signalfile = new File (signalfilename[f]);
			File tissuesfile = new File
				(Files.fileNamePrepend (signalfilename[f], "tissues_"));

			// Set up data set reader and writer.
			SignalDataSetReader reader =
				new SignalDataSetReader (signalfile);
			int H = reader.getHeight();
			int W = reader.getWidth();
			TissuesDataSetWriter writer =
				new TissuesDataSetWriter (tissuesfile, H, W);

			// Get time series.
			Series t_series = reader.getTimeSeries();
			int M = t_series.length();

			// Compute design matrix.
			double[][] A = new double [M] [N+1];
			for (int i = 0; i < M; ++ i)
				{
				double[] A_i = A[i];
				double t_i = t_series.x(i);
				for (int j = 0; j <= N; ++ j)
					{
					A_i[j] = SpinSignal.S (R1[j], t_i);
					}
				}

			// Analyze all pixels.
			int P = reader.getPixelCount();
			for (int i = 0; i < P; ++ i)
				{
				PixelSignal signal = reader.getPixelSignal (i);
				if (signal != null)
					{
					// Do the spin relaxometry analysis.
					PixelAnalysis.analyze
						(t_series, signal.S_measured(), R1_series, A,
						 rho_list, R1_list);

					// Write results to data set.
					writer.addPixelTissues
						(new PixelTissues (f, i, rho_list, R1_list));
					}
				}

			// All done.
			reader.close();
			writer.close();
			}

		// Stop timing.
		long t2 = System.currentTimeMillis();
		System.out.println ((t2-t1)+" msec");
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mri.SpinRelaxometrySeq <R1_lower> <R1_upper> <N> <signalfile> [<signalfile> ...]");
		System.err.println ("<R1_lower> = Lower R1 spin-lattice relaxation rate (1/sec)");
		System.err.println ("<R1_upper> = Upper R1 spin-lattice relaxation rate (1/sec)");
		System.err.println ("<N> = Number of R1 intervals");
		System.err.println ("<signalfile> = Input spin signal data set file");
		System.exit (1);
		}

	}
