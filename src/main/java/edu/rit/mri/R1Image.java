//******************************************************************************
//
// File:    R1Image.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.R1Image
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
 * Class R1Image is a program that displays an image depicting the locations of
 * tissues with certain spin-lattice relaxation rates in a magnetic resonance
 * image data set. A range of relaxation rates is specified on the command line.
 * The program displays an image. Pixels with no data are black. Pixels with no
 * tissues with relaxation rates in the given range are dark gray. Pixels with
 * one or more tissues with relaxation rates in the given range are medium gray
 * to white, depending on the spin densities; medium gray is a low spin density,
 * white is a high spin density. The program uses class {@linkplain
 * SignalDataSetReader} to read the data set.
 * <P>
 * Usage: java edu.rit.mri.test.R1Image <I>tissuesfile</I> <I>R1_lower</I>
 * <I>R1_upper</I>
 * <BR><I>tissuesfile</I> = Tissues data set file
 * <BR><I>R1_lower</I> = Lower <I>R</I>1 spin-lattice relaxation rate (1/sec)
 * <BR><I>R1_upper</I> = Upper <I>R</I>1 spin-lattice relaxation rate (1/sec)
 *
 * @author  Alan Kaminsky
 * @version 28-Jun-2008
 */
public class R1Image
	{

// Prevent construction.

	private R1Image()
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
		if (args.length != 3) usage();
		File tissuesfile = new File (args[0]);
		double R1_lower = Double.parseDouble (args[1]);
		double R1_upper = Double.parseDouble (args[2]);

		// Set up tissues data set reader.
		TissuesDataSetReader reader = new TissuesDataSetReader (tissuesfile);

		// Set up storage for image to be displayed.
		int H = reader.getHeight();
		int W = reader.getWidth();
		double[][] rhosum = new double [H] [W];
		byte[][] matrix = new byte [H] [W];
		PJGGrayImage image = new PJGGrayImage (H, W, matrix);

		// For each pixel, compute sum of spin densities of tissues having
		// R1_lower <= R1 < R1_upper.
		double rhomin = Double.POSITIVE_INFINITY;
		double rhomax = Double.NEGATIVE_INFINITY;
		for (int r = 0; r < H; ++ r)
			{
			for (int c = 0; c < W; ++ c)
				{
				PixelTissues tissues =
					reader.getPixelTissues (reader.indexFor (r, c));
				if (tissues != null)
					{
					image.setPixel (r, c, 0.25f);
					double rho = 0.0;
					for (int i = 0; i < tissues.numTissues(); ++ i)
						{
						double R1 = tissues.R1(i);
						if (R1_lower <= R1 && R1 <= R1_upper)
							{
							rho += tissues.rho(i);
							}
						}
					if (rho > 0.0)
						{
						rhomin = Math.min (rhomin, rho);
						rhomax = Math.max (rhomax, rho);
						rhosum[r][c] = rho;
						}
					}
				}
			}

		// Set the gray shade of each pixel.
		for (int r = 0; r < H; ++ r)
			{
			for (int c = 0; c < W; ++ c)
				{
				if (rhosum[r][c] > 0.0)
					{
					image.setPixel (r, c,
						(float)((rhosum[r][c]-rhomin)/(rhomax-rhomin)*0.5+0.5));
					}
				}
			}

		// Close tissues data set.
		reader.close();

		// Display plot.
		new DisplayableFrame
			("Tissues with "+args[1]+" <= R1 < "+args[2],
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
		System.err.println ("Usage: java edu.rit.mri.test.R1Image <tissuesfile> <R1_lower> <R1_upper>");
		System.err.println ("<tissuesfile> = Tissues data set file");
		System.err.println ("<R1_lower> = Lower R1 spin-lattice relaxation rate (1/sec)");
		System.err.println ("<R1_upper> = Upper R1 spin-lattice relaxation rate (1/sec)");
		System.exit (1);
		}

	}
