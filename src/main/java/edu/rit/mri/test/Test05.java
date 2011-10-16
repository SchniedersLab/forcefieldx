//******************************************************************************
//
// File:    Test05.java
// Package: edu.rit.mri.test
// Unit:    Class edu.rit.mri.test.Test05
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

import edu.rit.mri.PixelSignal;
import edu.rit.mri.SignalDataSetReader;
import edu.rit.mri.SpinSignal;
import edu.rit.mri.SpinSignalDifference;

import edu.rit.numeric.AggregateXYSeries;
import edu.rit.numeric.ListXYSeries;
import edu.rit.numeric.NonLinearLeastSquares;
import edu.rit.numeric.NonNegativeLeastSquares;
import edu.rit.numeric.Series;
import edu.rit.numeric.TooManyIterationsException;
import edu.rit.numeric.XYSeries;

import edu.rit.numeric.plot.Dots;
import edu.rit.numeric.plot.Plot;
import edu.rit.numeric.plot.Strokes;

import java.awt.Color;

import java.io.File;

import java.util.ArrayList;
import java.util.Collections;

/**
 * Class Test05 is a unit test main program for class {@linkplain
 * edu.rit.mri.SignalDataSetReader edu.rit.mri.SignalDataSetReader}. The program
 * reads the data for a specified pixel index. The program does a spin
 * relaxometry analysis on the data set. The program displays a plot of the data
 * set along with the reconstructed spin signal.
 * <P>
 * Usage: java edu.rit.mri.test.Test05 <I>signalfile</I> <I>R1_lower</I>
 * <I>R1_upper</I> <I>N</I> <I>index</I>
 * <BR><I>signalfile</I> = Input spin signal data set file
 * <BR><I>R1_lower</I> = Lower <I>R</I>1 spin-lattice relaxation rate (1/sec)
 * <BR><I>R1_upper</I> = Upper <I>R</I>1 spin-lattice relaxation rate (1/sec)
 * <BR><I>N</I> = Number of <I>R</I>1 intervals
 * <BR><I>index</I> = Pixel index
 *
 * @author  Alan Kaminsky
 * @version 24-Jun-2008
 */
public class Test05
	{

// Prevent construction.

	private Test05()
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
		Series tseries = reader.getTimeSeries();
		int M = tseries.length();
		PixelSignal signal = reader.getPixelSignal (index);
		if (signal == null)
			{
			System.err.println ("Test05: Pixel "+index+" has no data");
			System.exit (0);
			}
		Series dseries = signal.S_measured();
		AggregateXYSeries dataseries = new AggregateXYSeries (tseries, dseries);

//		// Print data X-Y series.
//		System.out.println ("Measured t, S(t)");
//		for (int i = 0; i < M; ++ i)
//			{
//			System.out.print (tseries.x(i));
//			System.out.print ('\t');
//			System.out.print (dseries.x(i));
//			System.out.println();
//			}

		// Do a spin relaxometry analysis using nonnegative linear least
		// squares. N spin relaxation rate intervals logarithmically spaced from
		// R1_lower to R1_upper.

		// Compute spin relaxation rates.
		double[] Rseries = new double [N+1];
		double logRlower = Math.log (R1_lower);
		double logRupper = Math.log (R1_upper);
		double interval = (logRupper - logRlower)/N;
		for (int j = 0; j <= N; ++ j)
			{
			Rseries[j] = Math.exp (logRlower + j*interval);
			}

		// Create nonnegative linear least squares solver.
		NonNegativeLeastSquares linsolver =
			new NonNegativeLeastSquares (M, N+1);

		// Find the solution.
		for (int i = 0; i < M; ++ i)
			{
			double[] a_i = linsolver.a[i];
			double t_i = dataseries.x(i);
			for (int j = 0; j <= N; ++ j)
				{
				a_i[j] = SpinSignal.S (Rseries[j], t_i);
				}
			linsolver.b[i] = dataseries.y(i);
			}
		linsolver.solve();
		double[] rhoseries = linsolver.x;

		// Print the solution.
		System.out.println ("Linear analysis: density, R1, T1");
		for (int j = 0; j <= N; ++ j)
			{
			if (rhoseries[j] > 0.0)
				{
				System.out.print (j);
				System.out.print (". ");
				System.out.print (rhoseries[j]);
				System.out.print ('\t');
				System.out.print (Rseries[j]);
				System.out.print ('\t');
				System.out.print (1.0/Rseries[j]);
				System.out.println();
				}
			}
		System.out.print ("chi^2 = ");
		System.out.print (linsolver.normsqr);
		System.out.println();

		// Find peaks in the solution. A peak occurs at index i if
		// rho[i] > rho[i-1] and rho[i] > rho[i+1].
		ArrayList<Double> rho = new ArrayList<Double>();
		ArrayList<Double> R = new ArrayList<Double>();
		ArrayList<Double> newrho = new ArrayList<Double>();
		ArrayList<Double> newR = new ArrayList<Double>();
		System.out.println ("Peaks in linear analysis: density, R1, T1");
		for (int j = 0; j <= N; ++ j)
			{
			if (rhoseries[j] > (j == 0 ? 0.0 : rhoseries[j-1]) &&
					rhoseries[j] > (j == N ? 0.0 : rhoseries[j+1]))
				{
				System.out.print (j);
				System.out.print (". ");
				System.out.print (rhoseries[j]);
				System.out.print ('\t');
				System.out.print (Rseries[j]);
				System.out.print ('\t');
				System.out.print (1.0/Rseries[j]);
				System.out.println();
				rho.add (rhoseries[j]);
				R.add (Rseries[j]);
				}
			}

		// Do a spin relaxometry analysis using nonlinear least squares. Peaks
		// in the linear analysis give the initial vector of densities and
		// rates.

		// Repeat until the solution is plausible.
		boolean plausible = false;
		int L = rho.size();
		while (L > 0 && ! plausible)
			{
			// Print solution before.
			System.out.println ("Nonlinear analysis, before: density, R1, T1");
			printSolution (dataseries, rho, R);
			
			// Create spin signal difference function. L = number of tissues.
			SpinSignalDifference fcn =
				new SpinSignalDifference (dataseries, L);

			// Create nonlinear least squares solver.
			NonLinearLeastSquares nonlinsolver =
				new NonLinearLeastSquares (fcn);

			// Find the solution.
			for (int i = 0; i < L; ++ i)
				{
				nonlinsolver.x[(i<<1)] = rho.get(i);
				nonlinsolver.x[(i<<1)+1] = R.get(i);
				}
			try
				{
				nonlinsolver.solve();
				newrho.clear();
				newR.clear();
				for (int i = 0; i < L; ++ i)
					{
					newrho.add (nonlinsolver.x[(i<<1)]);
					newR.add (nonlinsolver.x[(i<<1)+1]);
					}

				// Print solution after.
				System.out.println ("Nonlinear analysis, after: R1, T1");
				printSolution (dataseries, newrho, newR);

				// Decide if solution is plausible.
				plausible = checkPlausibility (dataseries, newrho, newR);
				}

			// Couldn't find a solution.
			catch (TooManyIterationsException exc)
				{
				System.out.println ("*** Too many iterations in nonlinear solver");
				plausible = false;
				}

			// If solution is not plausible, eliminate tissue with smallest
			// density and try again.
			if (! plausible)
				{
				double minrho = Double.MAX_VALUE;
				int mini = 0;
				for (int i = 0; i < L; ++ i)
					{
					if (rho.get(i) < minrho)
						{
						minrho = rho.get(i);
						mini = i;
						}
					}
				rho.remove (mini);
				R.remove (mini);
				L = rho.size();
				}

			// If solution is plausible, stop.
			else
				{
				rho = newrho;
				R = newR;
				}
			}

		// Create plot with data X-Y series.
		Plot plot = new Plot();
		plot.xAxisLength (1200)
			.xAxisMajorDivisions (20)
			.yAxisLength (600)
			.yAxisMajorDivisions (10)
			.seriesDots (Dots.circle (Color.BLACK, null, null, 5))
			.seriesColor (Color.BLACK)
			.seriesStroke (Strokes.solid (2))
			.xySeries (dataseries);

		// If we found a plausible solution, generate model spin signal X-Y
		// series and add to plot.
		if (L > 0)
			{
			ListXYSeries signalseries = new ListXYSeries();
			for (int i = 0; i < M; ++ i)
				{
				double t_i = dataseries.x(i);
				double s_i = 0.0;
				for (int j = 0; j < rho.size(); ++ j)
					{
					s_i += SpinSignal.S (rho.get(j), R.get(j), t_i);
					}
				signalseries.add (t_i, s_i);
				}
			plot.seriesDots (null)
				.seriesColor (Color.RED)
				.seriesStroke (Strokes.solid (1))
				.xySeries (signalseries);
			}
		else
			{
			System.out.println ("*** Could not find plausible solution");
			}

		// Display plot.
		plot.getFrame().setVisible (true);
		}

// Hidden operations.

	/**
	 * Print the given parameters.
	 */
	private static void printSolution
		(XYSeries dataseries,
		 ArrayList<Double> rho,
		 ArrayList<Double> R)
		{
		int L = rho.size();
		for (int i = 0; i < L; ++ i)
			{
			System.out.print (rho.get(i));
			System.out.print ('\t');
			System.out.print (R.get(i));
			System.out.print ('\t');
			System.out.print (1.0/R.get(i));
			System.out.println();
			}
		System.out.print ("chi^2 = ");
		System.out.print (computeChiSqr (dataseries, rho, R));
		System.out.println();
		}

	/**
	 * Calculate chi^2 between the measured data series and the model function
	 * determined by the given parameters.
	 */
	private static double computeChiSqr
		(XYSeries dataseries,
		 ArrayList<Double> rho,
		 ArrayList<Double> R)
		{
		int M = dataseries.length();
		int L = rho.size();
		double chisqr = 0.0;
		for (int i = 0; i < M; ++ i)
			{
			double t_i = dataseries.x(i);
			double s_i = 0.0;
			for (int j = 0; j < L; ++ j)
				{
				s_i += SpinSignal.S (rho.get(j), R.get(j), t_i);
				}
			double d = s_i - dataseries.y(i);
			chisqr += d*d;
			}
		return chisqr;
		}

	/**
	 * Decide if the given solution is plausible.
	 */
	private static boolean checkPlausibility
		(XYSeries dataseries,
		 ArrayList<Double> rho,
		 ArrayList<Double> R)
		{
		int M = dataseries.length();
		int L = rho.size();

		// If any density or rate is negative, solution is not plausible.
		for (int i = 0; i < L; ++ i)
			{
			if (rho.get(i) < 0.0)
				{
				System.out.println ("*** Negative density: "+rho.get(i));
				return false;
				}
			if (R.get(i) < 0.0)
				{
				System.out.println ("*** Negative rate: "+R.get(i));
				return false;
				}
			}

		// If relative difference between any two rates is too small, solution
		// is not plausible.
		for (int i = 0; i < L-1; ++ i)
			{
			double R_i = R.get(i);
			for (int j = i+1; j < L; ++ j)
				{
				double R_j = R.get(j);
				double reldiff = 2.0*Math.abs(R_i-R_j)/Math.abs(R_i+R_j);
				if (reldiff <= 0.001)
					{
					System.out.println
						("*** Rates too close: "+R_i+" and "+R_j+
						 ", reldiff = "+reldiff);
					return false;
					}
				}
			}

		// If sum of densities is too far from asymptotic measurement for large
		// t, solution is not plausible.
		double sumrho = 0.0;
		for (int i = 0; i < L; ++ i)
			{
			sumrho += rho.get(i);
			}
		double S_last = 0.0;
		int n = 0;
		for (int i = M-1; i >=0 && n < 7; -- i)
			{
			S_last += dataseries.y(i);
			++ n;
			}
		S_last /= n;
		double reldiff = Math.abs(sumrho-S_last)/Math.abs(S_last);
		if (reldiff >= 0.2)
			{
			System.out.println
				("*** Sum of densities = "+sumrho+
				 " doesn't agree with asymptotic measurement = "+S_last+
				 ", reldiff = "+reldiff);
			return false;
			}

		// Solution is plausible.
		return true;
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mri.test.Test05 <signalfile> <R1_lower> <R1_upper> <N> <index>");
		System.err.println ("<signalfile> = Input spin signal data set file");
		System.err.println ("<R1_lower> = Lower R1 spin-lattice relaxation rate (1/sec)");
		System.err.println ("<R1_upper> = Upper R1 spin-lattice relaxation rate (1/sec)");
		System.err.println ("<N> = Number of R1 intervals");
		System.err.println ("<index> = Pixel index");
		System.exit (1);
		}

	}
