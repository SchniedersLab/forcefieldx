//******************************************************************************
//
// File:    Test06.java
// Package: edu.rit.mri.test
// Unit:    Class edu.rit.mri.test.Test06
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

package edu.rit.mri.test;

import edu.rit.mri.PixelAnalysis;
import edu.rit.mri.PixelSignal;
import edu.rit.mri.SignalDataSetReader;
import edu.rit.mri.SpinSignal;

import edu.rit.numeric.AggregateXYSeries;
import edu.rit.numeric.ArraySeries;
import edu.rit.numeric.ListXYSeries;
import edu.rit.numeric.Series;

import edu.rit.numeric.plot.Dots;
import edu.rit.numeric.plot.Plot;
import edu.rit.numeric.plot.Strokes;

import java.awt.Color;

import java.io.File;

import java.util.ArrayList;
import java.util.List;

/**
 * Class Test06 is a unit test main program for class {@linkplain
 * edu.rit.mri.PixelAnalysis edu.rit.mri.PixelAnalysis}. Using an instance of
 * class {@linkplain edu.rit.mri.SignalDataSetReader
 * edu.rit.mri.SignalDataSetReader}, the program reads the data for a specified
 * pixel index. The program does a spin relaxometry analysis on the data set.
 * The program displays a plot of the data set along with the reconstructed spin
 * signal.
 * <P>
 * Usage: java edu.rit.mri.test.Test06 <I>signalfile</I> <I>R1_lower</I>
 * <I>R1_upper</I> <I>N</I> <I>index</I>
 * <BR><I>signalfile</I> = Input spin signal data set file
 * <BR><I>R1_lower</I> = Lower <I>R</I>1 spin-lattice relaxation rate (1/sec)
 * <BR><I>R1_upper</I> = Upper <I>R</I>1 spin-lattice relaxation rate (1/sec)
 * <BR><I>N</I> = Number of <I>R</I>1 intervals
 * <BR><I>index</I> = Pixel index
 *
 * @author  Alan Kaminsky
 * @version 21-Jun-2008
 */
public class Test06
	{

// Prevent construction.

	private Test06()
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
		if (args.length != 5) usage();
		File signalfile = new File (args[0]);
		double R1_lower = Double.parseDouble (args[1]);
		double R1_upper = Double.parseDouble (args[2]);
		int N = Integer.parseInt (args[3]);
		int index = Integer.parseInt (args[4]);

		// Set up spin signal data set reader.
		SignalDataSetReader reader = new SignalDataSetReader (signalfile);

		// Get data X-Y series.
		Series t_series = reader.getTimeSeries();
		int M = t_series.length();
		PixelSignal signal = reader.getPixelSignal (index);
		if (signal == null)
			{
			System.err.println ("Test06: Pixel "+index+" has no data");
			System.exit (0);
			}
		Series S_series = signal.S_measured();

//		// Print data X-Y series.
//		System.out.println ("Measured t, S(t)");
//		for (int i = 0; i < M; ++ i)
//			{
//			System.out.print (tseries.x(i));
//			System.out.print ('\t');
//			System.out.print (dseries.x(i));
//			System.out.println();
//			}

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

		// Do the spin relaxometry analysis.
		ArrayList<Double> rho_list = new ArrayList<Double>();
		ArrayList<Double> R1_list = new ArrayList<Double>();
		PixelAnalysis.analyze
			(t_series, S_series, R1_series, A, rho_list, R1_list);

		// Print solution.
		printSolution (t_series, S_series, rho_list, R1_list);

		// Create plot with data X-Y series.
		Plot plot = new Plot();
		plot.xAxisLength (1200)
			.xAxisMajorDivisions (20)
			.yAxisLength (600)
			.yAxisMajorDivisions (10)
			.seriesDots (Dots.circle (Color.BLACK, null, null, 5))
			.seriesColor (Color.BLACK)
			.seriesStroke (Strokes.solid (2))
			.xySeries (new AggregateXYSeries (t_series, S_series));

		// If we found a plausible solution, generate model spin signal X-Y
		// series and add to plot.
		int L = rho_list.size();
		if (L > 0)
			{
			ListXYSeries signal_series = new ListXYSeries();
			for (int i = 0; i < M; ++ i)
				{
				double t_i = t_series.x(i);
				double s_i = 0.0;
				for (int j = 0; j < L; ++ j)
					{
					s_i += SpinSignal.S (rho_list.get(j), R1_list.get(j), t_i);
					}
				signal_series.add (t_i, s_i);
				}
			plot.seriesDots (null)
				.seriesColor (Color.RED)
				.seriesStroke (Strokes.solid (1))
				.xySeries (signal_series);
			}

		// Display plot.
		plot.getFrame().setVisible (true);
		}

// Hidden operations.

	/**
	 * Print the given solution.
	 */
	private static void printSolution
		(Series t_series,
		 Series S_series,
		 List<Double> rho_list,
		 List<Double> R1_list)
		{
		int M = t_series.length();
		int L = rho_list.size();
		if (L == 0)
			{
			System.out.println ("Test06: Could not find solution");
			}
		else
			{
			System.out.println ("rho\tR1\tT1");
			for (int i = 0; i < L; ++ i)
				{
				System.out.print (rho_list.get(i));
				System.out.print ('\t');
				System.out.print (R1_list.get(i));
				System.out.print ('\t');
				System.out.print (1.0/R1_list.get(i));
				System.out.println();
				}

			double chisqr =
				computeChiSqr (t_series, S_series, rho_list, R1_list);
			System.out.print ("chi^2 = ");
			System.out.print (chisqr);
			System.out.println();

			int dof = M - 2*L;
			System.out.print ("dof   = ");
			System.out.print (dof);
			System.out.println();

			double sigma = Math.sqrt(chisqr/dof);
			System.out.print ("sigma = ");
			System.out.print (sigma);
			System.out.println();
			}
		}

	/**
	 * Calculate chi^2 between the measured data series and the model function
	 * determined by the given parameters.
	 */
	private static double computeChiSqr
		(Series t_series,
		 Series S_series,
		 List<Double> rho_list,
		 List<Double> R1_list)
		{
		int M = t_series.length();
		int L = rho_list.size();
		double chisqr = 0.0;
		for (int i = 0; i < M; ++ i)
			{
			double t_i = t_series.x(i);
			double s_i = 0.0;
			for (int j = 0; j < L; ++ j)
				{
				s_i += SpinSignal.S (rho_list.get(j), R1_list.get(j), t_i);
				}
			double d = s_i - S_series.x(i);
			chisqr += d*d;
			}
		return chisqr;
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mri.test.Test06 <signalfile> <R1_lower> <R1_upper> <N> <index>");
		System.err.println ("<signalfile> = Input spin signal data set file");
		System.err.println ("<R1_lower> = Lower R1 spin-lattice relaxation rate (1/sec)");
		System.err.println ("<R1_upper> = Upper R1 spin-lattice relaxation rate (1/sec)");
		System.err.println ("<N> = Number of R1 intervals");
		System.err.println ("<index> = Pixel index");
		System.exit (1);
		}

	}
