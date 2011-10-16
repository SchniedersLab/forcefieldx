//******************************************************************************
//
// File:    SpinSignalImage.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.SpinSignalImage
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

import edu.rit.image.PJGGrayImage;

import edu.rit.numeric.AggregateXYSeries;
import edu.rit.numeric.Interpolation;
import edu.rit.numeric.Series;

import edu.rit.swing.DisplayableFrame;

import java.io.File;

import javax.swing.JFrame;

/**
 * Class SpinSignalImage is a program that displays an image depicting the spin
 * signals at a given time value for a magnetic resonance image data set. Each
 * pixel's gray shade represents the value of the pixel's spin signal, with
 * &minus;32768 being black and +32768 being white. The program uses class
 * {@linkplain SignalDataSetReader} to read the data set.
 * <P>
 * Usage: java edu.rit.mri.SpinSignalImage <I>signalfile</I> <I>t</I>
 * <BR><I>signalfile</I> = Spin signal data set file
 * <BR><I>t</I> = Time (sec)
 *
 * @author  Alan Kaminsky
 * @version 21-Jun-2008
 */
public class SpinSignalImage
	{

// Prevent construction.

	private SpinSignalImage()
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
		if (args.length != 2) usage();
		File signalfile = new File (args[0]);
		double t = Double.parseDouble (args[1]);

		// Set up spin signal data set reader.
		SignalDataSetReader reader = new SignalDataSetReader (signalfile);

		// Set up image to be displayed.
		int H = reader.getHeight();
		int W = reader.getWidth();
		byte[][] matrix = new byte [H] [W];
		PJGGrayImage image = new PJGGrayImage (H, W, matrix);

		// Get time series.
		Series t_series = reader.getTimeSeries();

		// Set the gray shade of each pixel.
		for (int r = 0; r < H; ++ r)
			{
			for (int c = 0; c < W; ++ c)
				{
				// Get spin signal series.
				PixelSignal signal =
					reader.getPixelSignal (reader.indexFor (r, c));
				if (signal != null)
					{
					Series S_series = signal.S_measured();

					// Interpolate spin signal at the given time value.
					double S =
						new Interpolation
							(new AggregateXYSeries (t_series, S_series))
							.f(t);

					// Scale range -32768.0 .. +32768.0 to range 0.0 .. 1.0 and
					// set gray shade.
					image.setPixel (r, c, (float)((S + 32768.0)/65536.0));
					}
				}
			}

		// Close spin signal data set.
		reader.close();

		// Display plot.
		new DisplayableFrame
			("Spin Signal, t = "+t,
			 image.getDisplayable(),
			 JFrame.EXIT_ON_CLOSE)
			.setVisible (true);
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mri.SpinSignalImage <signalfile> <t>");
		System.err.println ("<signalfile> = Spin signal data set file");
		System.err.println ("<t> = Time (sec)");
		System.exit (1);
		}

	}
