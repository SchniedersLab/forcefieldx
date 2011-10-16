//******************************************************************************
//
// File:    R1Distribution.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.R1Distribution
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

import edu.rit.numeric.ListXYSeries;

import edu.rit.numeric.plot.Plot;

import edu.rit.numeric.plot.impl.LinearAxis;

import java.io.File;

import java.util.ArrayList;
import java.util.Arrays;

import java.text.DecimalFormat;

/**
 * Class R1Distribution is a main program that displays the results of a spin
 * relaxometry analysis of a magnetic resonance image.
 * <P>
 * The program displays a plot of the cumulative distribution of the computed
 * spin-lattice relaxation rates <I>R</I>1 for all the tissues. The program
 * reads the relaxation rates from a tissues data set file using class
 * {@linkplain TissuesDataSetReader}. The program also prints the deciles of the
 * distribution. The lower and upper bounds of the relaxation rates to plot may
 * be specified on the command line; if omitted, the default is to include all
 * the relaxation rates in the data set.
 * <P>
 * Usage: java edu.rit.mri.R1Distribution <I>tissuesfile</I> [ <I>R1_lower</I>
 * <I>R1_upper</I> ]
 * <BR><I>tissuesfile</I> = Tissues data set file
 * <BR><I>R1_lower</I> = Lower <I>R</I>1 spin-lattice relaxation rate (1/sec)
 * <BR><I>R1_upper</I> = Upper <I>R</I>1 spin-lattice relaxation rate (1/sec)
 *
 * @author  Alan Kaminsky
 * @version 21-Jun-2008
 */
public class R1Distribution
	{

// Prevent construction.

	private R1Distribution()
		{
		}

// Main program.

	static final DecimalFormat FMT1 = new DecimalFormat ("0.0");
	static final DecimalFormat FMT2 = new DecimalFormat ("0.0E0");

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Parse command line arguments.
		if (args.length != 1 && args.length != 3) usage();
		File tissuesfile = new File (args[0]);
		double R1_lower = -1.0;
		double R1_upper = -1.0;
		if (args.length == 3)
			{
			R1_lower = Double.parseDouble (args[1]);
			R1_upper = Double.parseDouble (args[2]);
			}

		// Make a list of all spin-lattice relaxation rates in the data set.
		TissuesDataSetReader reader = new TissuesDataSetReader (tissuesfile);
		ArrayList<Double> R1_list = new ArrayList<Double>();
		int P = reader.getPixelCount();
		for (int i = 0; i < P; ++ i)
			{
			PixelTissues tissues = reader.getPixelTissues(i);
			if (tissues != null)
				{
				int L = tissues.numTissues();
				for (int j = 0; j < L; ++ j)
					{
					R1_list.add (tissues.R1(j));
					}
				}
			}

		// Convert the list of relaxation rates to a sorted array.
		int N = R1_list.size();
		double[] R1_array = new double [N];
		for (int i = 0; i < N; ++ i)
			{
			R1_array[i] = R1_list.get(i);
			}
		Arrays.sort (R1_array);

		// Print deciles.
		System.out.println ("Decile\tR1");
		for (int i = 0; i <= 9; ++ i)
			{
			System.out.print (FMT1.format (i/10.0));
			System.out.print ('\t');
			System.out.print (R1_array[N*i/10]);
			System.out.println();
			}
		System.out.print ("1.0\t");
		System.out.print (R1_array[N-1]);
		System.out.println();

		// Make an X-Y series of the cumulative distribution.
		ListXYSeries series = new ListXYSeries();
		series.add (0, 0);
		for (int i = 0; i < N; ++ i)
			{
			series.add (R1_array[i], ((double) i)/((double) N));
			series.add (R1_array[i], ((double)(i+1))/((double) N));
			}

		// Create plot of cumulative distribution.
		if (R1_lower == -1.0)
			{
			R1_lower = 0.0;
			R1_upper = LinearAxis.autoscale (R1_array[N-1]);
			}
		Plot plot = new Plot();
		plot.plotTitle ("<I>R</I>1 Distribution")
			.xAxisLength (800)
			.xAxisStart (R1_lower)
			.xAxisEnd (R1_upper)
			.xAxisMajorDivisions (20)
			.xAxisTickFormat (FMT2)
			.xAxisTitle ("Spin-Lattice Relaxation Rate, <I>R</I>1 (1/sec)")
			.yAxisLength (400)
			.yAxisStart (0)
			.yAxisEnd (1)
			.yAxisMajorDivisions (10)
			.yAxisTickFormat (FMT1)
			.yAxisTitle ("Cumulative Distribution")
			.seriesDots (null)
			.xySeries (series);

		// Display plot.
		plot.getFrame().setVisible (true);
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mri.R1Distribution <tissuesfile> [<R1_lower> <R1_upper>]");
		System.err.println ("<tissuesfile> = Tissues data set file");
		System.err.println ("<R1_lower> = Lower R1 spin-lattice relaxation rate (1/sec)");
		System.err.println ("<R1_upper> = Upper R1 spin-lattice relaxation rate (1/sec)");
		System.exit (1);
		}

	}
