//******************************************************************************
//
// File:    PlotPixel.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.PlotPixel
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

import edu.rit.numeric.AggregateXYSeries;
import edu.rit.numeric.Series;

import edu.rit.numeric.plot.Dots;
import edu.rit.numeric.plot.Plot;
import edu.rit.numeric.plot.Strokes;

import edu.rit.numeric.plot.impl.LinearAxis;

import java.awt.Color;

import java.io.File;

import java.text.DecimalFormat;

/**
 * Class PlotPixel is a main program that displays the results of a spin
 * relaxometry analysis of a magnetic resonance image.
 * <P>
 * The program displays a plot of the measured spin signal for the pixel at a
 * given index. The program reads the spin signal from a spin signal data set
 * file using class {@linkplain SignalDataSetReader}.
 * <P>
 * The program includes a plot of the fitted spin signal for the pixel as
 * determined by the spin relaxometry analysis. The program reads the analysis
 * results from a tissues data set file using class {@linkplain
 * TissuesDataSetReader}. The program also prints the tissues' computed spin
 * densities <I>&rho;</I><SUB><I>j</I></SUB>, spin-lattice relaxation rates
 * <I>R</I>1<SUB><I>j</I></SUB>, and spin-lattice relaxation times
 * <I>T</I>1<SUB><I>j</I></SUB>. The program prints <I>&chi;</I><SUP>2</SUP>,
 * the sum of the squared differences between the measured spin signal and the
 * fitted spin signal; the degrees of freedom (dof) in <I>&chi;</I><SUP>2</SUP>;
 * and the standard deviation <I>&sigma;</I> of the errors in the measured
 * signal derived from <I>&chi;</I><SUP>2</SUP>.
 * <P>
 * If no tissues data set file is specified, the analysis results are not
 * plotted or printed.
 * <P>
 * Usage: java edu.rit.mri.PlotPixel <I>signalfile</I> [ <I>tissuesfile</I> ]
 * <I>index</I>
 * <BR><I>signalfile</I> = Spin signal data set file
 * <BR><I>tissuesfile</I> = Tissues data set file
 * <BR><I>index</I> = Pixel index
 *
 * @author  Alan Kaminsky
 * @version 21-Jun-2008
 */
public class PlotPixel
	{

// Prevent construction.

	private PlotPixel()
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
		File signalfile = null;
		File tissuesfile = null;
		int index = 0;
		if (args.length == 2)
			{
			signalfile = new File (args[0]);
			index = Integer.parseInt (args[1]);
			}
		else if (args.length == 3)
			{
			signalfile = new File (args[0]);
			tissuesfile = new File (args[1]);
			index = Integer.parseInt (args[2]);
			}
		else
			{
			usage();
			}

		// Read time series and pixel signal for the given pixel.
		SignalDataSetReader signalReader = new SignalDataSetReader (signalfile);
		Series t_series = signalReader.getTimeSeries();
		PixelSignal signal = signalReader.getPixelSignal (index);
		signalReader.close();
		if (signal == null)
			{
			System.err.println ("PlotPixel: Pixel "+index+" has no data");
			System.exit (0);
			}
		Series S_measured = signal.S_measured();
		int M = t_series.length();

		// Read pixel tissues for the given pixel.
		int L = 0;
		Series S_computed = null;
		if (tissuesfile != null)
			{
			TissuesDataSetReader tissuesReader =
				new TissuesDataSetReader (tissuesfile);
			PixelTissues tissues = tissuesReader.getPixelTissues (index);
			tissuesReader.close();
			if (tissues == null)
				{
				System.out.println ("No solution found");
				}
			else
				{
				L = tissues.numTissues();
				if (L == 0)
					{
					System.out.println ("No solution found");
					}
				else
					{
					// Get computed spin signal.
					S_computed = tissues.S_series (t_series);

					// Print analysis results.
					System.out.println ("rho\tR1\tT1");
					for (int i = 0; i < L; ++ i)
						{
						System.out.print (tissues.rho(i));
						System.out.print ('\t');
						System.out.print (tissues.R1(i));
						System.out.print ('\t');
						System.out.print (1.0/tissues.R1(i));
						System.out.println();
						}

					double chisqr = 0.0;
					for (int i = 0; i < M; ++ i)
						{
						double d = S_measured.x(i) - S_computed.x(i);
						chisqr += d*d;
						}
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
			}

		// Create plot with measured spin signal X-Y series.
		double maxY = 0.0;
		maxY = Math.max (maxY, Math.abs (S_measured.minX()));
		maxY = Math.max (maxY, Math.abs (S_measured.maxX()));
		if (L > 0)
			{
			maxY = Math.max (maxY, Math.abs (S_computed.minX()));
			maxY = Math.max (maxY, Math.abs (S_computed.maxX()));
			}
		maxY = LinearAxis.autoscale (maxY);
		Plot plot = new Plot();
		plot.plotTitle ("Spin Signal vs. Time for Pixel "+index)
			.leftMargin (72)
			.xAxisLength (576)
			.xAxisMajorDivisions (20)
			.yAxisLength (288)
			.xAxisTitle ("Time, <I>t</I> (sec)")
			.yAxisMajorDivisions (10)
			.yAxisStart (-maxY)
			.yAxisEnd (+maxY)
			.yAxisTickFormat (new DecimalFormat ("0.0E0"))
			.yAxisTitle ("Spin signal, <I>S</I> (<I>t</I>)")
			.yAxisTitleOffset (54)
			.seriesDots (Dots.circle (Color.BLACK, null, null, 5))
			.seriesColor (Color.BLACK)
			//.seriesStroke (Strokes.solid (1))
			.seriesStroke (null)
			.xySeries (new AggregateXYSeries (t_series, S_measured));

		// If we found a solution, add computed spin signal X-Y series to plot.
		if (L > 0)
			{
			plot.seriesDots (null)
				.seriesColor (Color.BLACK)
				.seriesStroke (Strokes.solid (2))
				.xySeries (new AggregateXYSeries (t_series, S_computed));
			}

		// Display plot.
		plot.getFrame().setVisible (true);
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mri.PlotPixel <signalfile> [<tissuesfile>] <index>");
		System.err.println ("<signalfile> = Spin signal data set file");
		System.err.println ("<tissuesfile> = Tissues data set file");
		System.err.println ("<index> = Pixel index");
		System.exit (1);
		}

	}
