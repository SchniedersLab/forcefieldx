//******************************************************************************
//
// File:    R1Histogram.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.R1Histogram
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

import edu.rit.numeric.plot.Plot;

import java.io.File;

import java.text.DecimalFormat;

/**
 * Class R1Histogram is a main program that displays the results of a spin
 * relaxometry analysis of a magnetic resonance image.
 * <P>
 * The program displays a histogram of the computed spin-lattice relaxation
 * rates <I>R</I>1 for all the tissues. The program reads the relaxation rates
 * from a tissues data set file using class {@linkplain TissuesDataSetReader}.
 * The relaxation rate lower and upper bounds and the number of intervals for
 * the histogram are specified on the command line.
 * <P>
 * Usage: java edu.rit.mri.R1Histogram <I>tissuesfile</I> <I>R1_lower</I>
 * <I>R1_upper</I> <I>N</I>
 * <BR><I>tissuesfile</I> = Tissues data set file
 * <BR><I>R1_lower</I> = Lower <I>R</I>1 spin-lattice relaxation rate (1/sec)
 * <BR><I>R1_upper</I> = Upper <I>R</I>1 spin-lattice relaxation rate (1/sec)
 * <BR><I>N</I> = Number of intervals
 *
 * @author  Alan Kaminsky
 * @version 27-Jun-2008
 */
public class R1Histogram
	{

// Prevent construction.

	private R1Histogram()
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
		if (args.length != 4) usage();
		File tissuesfile = new File (args[0]);
		double R1_lower = Double.parseDouble (args[1]);
		double R1_upper = Double.parseDouble (args[2]);
		int N = Integer.parseInt (args[3]);

		double binWidth = (R1_upper - R1_lower)/N;

		// Set up histogram bins.
		double[] x = new double [N+1];
		double[] y = new double [N+1];
		for (int i = 0; i <= N; ++ i)
			{
			x[i] = R1_lower + i*binWidth;
			}

		// Count all spin-lattice relaxation rates in the data set.
		TissuesDataSetReader reader = new TissuesDataSetReader (tissuesfile);
		int P = reader.getPixelCount();
		for (int i = 0; i < P; ++ i)
			{
			PixelTissues tissues = reader.getPixelTissues(i);
			if (tissues != null)
				{
				int L = tissues.numTissues();
				for (int j = 0; j < L; ++ j)
					{
					double R1 = tissues.R1(j);
					int k = (int) Math.floor((R1 - R1_lower)/binWidth);
					if (0 <= k && k <= N) y[k] += 1.0;
					}
				}
			}

		// Create plot of histogram.
		Plot plot = new Plot();
		plot.plotTitle ("<I>R</I><SUB>1</SUB> Histogram")
			.leftMargin (54)
			.xAxisLength (800)
			.xAxisStart (R1_lower)
			.xAxisEnd (R1_upper)
			.xAxisMajorDivisions (20)
			.xAxisTickFormat (FMT2)
			.xAxisTitle ("Spin-lattice relaxation rate, <I>R</I><SUB>1</SUB> (sec<SUP>\u20121</SUP>)")
			.yAxisLength (400)
			.yAxisMajorDivisions (10)
			.yAxisTickFormat (FMT2)
			.yAxisTitle ("Count")
			.yAxisTitleOffset (42)
			.xySeries (x, y);

		// Display plot.
		plot.getFrame().setVisible (true);
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mri.R1Histogram <tissuesfile> <R1_lower> <R1_upper> <N>");
		System.err.println ("<tissuesfile> = Tissues data set file");
		System.err.println ("<R1_lower> = Lower R1 spin-lattice relaxation rate (1/sec)");
		System.err.println ("<R1_upper> = Upper R1 spin-lattice relaxation rate (1/sec)");
		System.err.println ("<N> = Number of intervals");
		System.exit (1);
		}

	}
