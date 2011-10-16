//******************************************************************************
//
// File:    RenderSeq.java
// Package: edu.rit.hyb.antimatter
// Unit:    Class edu.rit.hyb.antimatter.RenderSeq
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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

package edu.rit.hyb.antimatter;

import edu.rit.pj.Comm;

import edu.rit.util.Random;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;

import java.awt.geom.Ellipse2D;

import java.awt.image.BufferedImage;
import java.awt.image.IndexColorModel;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;

import javax.imageio.ImageIO;

/**
 * Class RenderSeq is a sequential program that renders a visualization of the
 * Antimatter Simulation. This program's chief purpose is to measure the
 * rendering performance. The program does the following:
 * <P>
 * Repeat for <I>F</I> frames:
 * <BR>&nbsp;&nbsp;&nbsp;&nbsp;Place <I>N</I> antiprotons at random positions in
 * the square from (0,0) to (<I>W,W</I>).
 * <BR>&nbsp;&nbsp;&nbsp;&nbsp;Create a <I>W</I>x<I>W</I>-pixel image of the
 * antiprotons.
 * <BR>&nbsp;&nbsp;&nbsp;&nbsp;Store the image in a PNG file.
 * <P>
 * Usage: java edu.rit.hyb.antimatter.RenderSeq <I>seed</I> <I>N</I> <I>W</I>
 * <I>F</I> <I>file</I>
 * <BR><I>seed</I> = Random seed
 * <BR><I>N</I> = Number of antiprotons
 * <BR><I>W</I> = Image is <I>W</I>x<I>W</I> pixels
 * <BR><I>F</I> = Number of frames (files)
 * <BR><I>file</I> = Files are named "<I>file</I>_0000.png",
 * "<I>file</I>_0001.png", etc.
 * <P>
 * Class RenderSeq executes sequentially in a single thread.
 *
 * @author  Alan Kaminsky
 * @version 31-Mar-2007
 */
public class RenderSeq
	{

// Prevent construction.

	private RenderSeq()
		{
		}

// Hidden constants.

	private static final double DIAM = 2.0;
	private static final double DIAM_OVER_2 = DIAM/2.0;
	//private static final Color BACKGROUND_COLOR = new Color (0.9f, 0.9f, 0.9f);
	private static final Color BACKGROUND_COLOR = Color.black;
	private static final Color ANTIPROTON_COLOR = Color.red;

// Hidden data members.

	// Command line arguments.
	static long seed;
	static int N;
	static int W;
	static int F;
	static String file;

	// Pseudorandom number generator.
	static Random prng;

	// For drawing antiprotons.
	static Ellipse2D dot = new Ellipse2D.Double();
	static IndexColorModel colormodel;
	static BufferedImage image;
	static Graphics2D g2d;

	// For generating file names.
	static StringBuilder filename = new StringBuilder();

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		Comm.init (args);

		// Start timing.
		long t1 = System.currentTimeMillis();

		// Parse command line arguments.
		if (args.length != 5) usage();
		seed = Long.parseLong (args[0]);
		N = Integer.parseInt (args[1]);
		W = Integer.parseInt (args[2]);
		F = Integer.parseInt (args[3]);
		file = args[4];

		// Set up PRNG.
		prng = Random.getInstance (seed);

		// Set up for drawing.
		byte[] red = new byte [256];
		byte[] green = new byte [256];
		byte[] blue = new byte [256];
		for (int i = 0; i < 256; ++ i)
			{
			red[i] = (byte) i;
			}
		colormodel = new IndexColorModel (8, 256, red, green, blue);
		image = new BufferedImage
			(W, W, BufferedImage.TYPE_BYTE_INDEXED, colormodel);
		g2d = image.createGraphics();
		g2d.setRenderingHint
			(RenderingHints.KEY_ANTIALIASING,
			 RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setRenderingHint
			(RenderingHints.KEY_DITHERING,
			 RenderingHints.VALUE_DITHER_DISABLE);

		// Create F frames.
		for (int f = 0; f < F; ++ f)
			{
			// Fill in the background.
			g2d.setColor (BACKGROUND_COLOR);
			g2d.fillRect (0, 0, W, W);

			// Draw N particles at random positions from (0,0) to (R,R).
			g2d.setColor (ANTIPROTON_COLOR);
			for (int i = 0; i < N; ++ i)
				{
				dot.setFrame
					(W * prng.nextDouble() - DIAM_OVER_2,
					 W * prng.nextDouble() - DIAM_OVER_2,
					 DIAM, DIAM);
				g2d.fill (dot);
				}

			// Generate file name.
			filename.setLength (0);
			filename.append (f);
			while (filename.length() < 4) filename.insert (0, '0');
			filename.insert (0, '_');
			filename.insert (0, file);
			filename.append (".png");

			// Write PNG file.
			ImageIO.write
				(image,
				 "png",
				 new BufferedOutputStream
					(new FileOutputStream
						(filename.toString())));
			}

		// Stop timing.
		long t2 = System.currentTimeMillis();
		System.out.println ((t2-t1) + " msec");
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.hyb.antimatter.RenderSeq <seed> <N> <W> <F> <file>");
		System.err.println ("<seed> = Random seed");
		System.err.println ("<N> = Number of antiprotons");
		System.err.println ("<W> = Image is <W>x<W> pixels");
		System.err.println ("<F> = Number of frames (files)");
		System.err.println ("<file> = Files are named \"<file>_0000.png\", \"<file>_0001.png\", etc.");
		System.exit (1);
		}

	}
