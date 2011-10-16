//******************************************************************************
//
// File:    AntiprotonHyb.java
// Package: edu.rit.hyb.antimatter
// Unit:    Class edu.rit.hyb.antimatter.AntiprotonHyb
//
// This Java source file is copyright (C) 2006 by Alan Kaminsky. All rights
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

import edu.rit.mp.DoubleBuf;

import edu.rit.pj.BarrierAction;
import edu.rit.pj.Comm;
import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import edu.rit.util.Random;
import edu.rit.util.Range;

import edu.rit.vector.Vector2D;

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
 * Class AntiprotonHyb is a hybrid SMP cluster parallel program that calculates
 * the positions of the antiprotons as a function of time.
 * <P>
 * Usage: java -Dpj.np=<I>Kp</I> -Dpj.nt=<I>Kt</I>
 * edu.rit.hyb.antimatter.AntiprotonHyb <I>seed</I> <I>N</I> <I>R</I>
 * <I>frames</I> <I>steps</I> <I>dt</I> <I>W</I> <I>file</I>
 * <P>
 * The program runs in <I>Kp</I> parallel processes. Within each process, the
 * program runs in <I>Kt</I> parallel threads. The program:
 * <OL TYPE=1>
 * <LI>
 * Initializes a pseudorandom number generator with <I>seed</I>.
 * <LI>
 * Generates <I>N</I> antiprotons positioned at random in the square from
 * (0,0) to (<I>R,R</I>).
 * <LI>
 * Sets each antiproton's initial velocity to 0.
 * <LI>
 * Calculates <I>frames</I> visualization frames. For each frame:
 * <OL TYPE=a>
 * <LI>
 * Performs <I>steps</I> time steps, calculating the antiprotons' positions
 * after an elapsed time of <I>dt</I>.
 * <LI>
 * Renders a <I>W</I>x<I>W</I>-pixel image of the antiprotons' positions.
 * <LI>
 * Stores the image in a PNG file.
 * </OL>
 * </OL>
 * The image files are named <TT>"<I>file</I>_0000.png"</TT>,
 * <TT>"<I>file</I>_0001.png"</TT>, and so on. The first file depicts the
 * antiprotons' initial positions. Each subsequent file depicts the antiprotons'
 * positions <I>steps</I> time steps after the previous file.
 * <P>
 * The computation is performed in parallel in multiple processes. One process
 * computes the visualization. The remaining <I>Kp</I>-1 processes compute the
 * antiproton positions in parallel, with <I>Kt</I> parallel threads in each
 * process. Thus, <I>Kp</I> must be 2 or greater. The program measures the
 * computation's running time.
 * <P>
 * Each process holds the entire position array, one slice of the velocity
 * array, and one slice of the force array. At each time step, each process
 * calculates its own slice of the force array using the entire position array.
 * Each process uses its slice of the force array to update its slices of the
 * velocity and position arrays. Then the processes do an all-gather of the
 * position array slices, so every processor has the complete new position array
 * for the next time step.
 *
 * @author  Alan Kaminsky
 * @version 11-Nov-2007
 */
public class AntiprotonHyb
	{

// Prevent construction.

	private AntiprotonHyb()
		{
		}

// Hidden constants.

	// Charge on an antiproton.
	static final double QP = 3.0;

	// Charge per unit length on the trap.
	static final double QT = 3.0;

	// For visualization.
	static final double DIAM = 2.0;
	static final double DIAM_OVER_2 = DIAM/2.0;
	static final Color BACKGROUND_COLOR = Color.black;
	static final Color ANTIPROTON_COLOR = Color.red;

// Program variables.

	// World communicator.
	static Comm world;
	static int size;
	static int rank;

	// Command line arguments.
	static long seed;
	static int N;
	static double R;
	static int frames;
	static int steps;
	static double dt;
	static int W;
	static String file;

	static double onehalfdtsqr;
	static double QP_QP;
	static double QP_QT;

	// Communicator for processes calculating antiproton positions.
	static Comm p_comm;
	static int p_size;
	static int p_rank;

	// Array slices.
	static Range[] ranges;
	static Range myrange;
	static int mylb;
	static int myub;
	static int mylength;

	// Force, velocity, and position (x,y) arrays.
	static Vector2D[] f;
	static Vector2D[] v;
	static Vector2D[] p;

	// Communication buffers.
	static DoubleBuf[] p_slices;
	static DoubleBuf myp_slice;

	// For drawing antiprotons.
	static double scale;
	static Ellipse2D dot;
	static IndexColorModel colormodel;
	static BufferedImage image;
	static Graphics2D g2d;

	// For generating file names.
	static StringBuilder filename = new StringBuilder();

	// For thread parallelism within each process.
	static ParallelTeam team;
	static ParallelRegion region;

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Start timing.
		long time = -System.currentTimeMillis();

		// Initialize middleware.
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();

		// Parse command line arguments.
		if (args.length != 8) usage();
		seed = Long.parseLong (args[0]);
		N = Integer.parseInt (args[1]);
		R = Double.parseDouble (args[2]);
		frames = Integer.parseInt (args[3]);
		steps = Integer.parseInt (args[4]);
		dt = Double.parseDouble (args[5]);
		W = Integer.parseInt (args[6]);
		file = args[7];

		onehalfdtsqr = 0.5 * dt * dt;
		QP_QP = QP * QP;
		QP_QT = QP * QT * N;

		// Set up position communicator consisting of processes 0 .. size-2.
		p_comm = world.createComm (rank <= size-2);

		// Processes 0 .. size-2 do the antiproton position computations,
		// process size-1 does the visualization computations.
		if (rank <= size-2)
			{
			computePositions();
			}
		else
			{
			computeVisualizations();
			}

		// Stop timing.
		time += System.currentTimeMillis();
		System.out.println (time + " msec " + rank);
		}

// Hidden operations.

	/**
	 * Antiproton position computations in processes 0 .. size-2.
	 */
	private static void computePositions()
		throws Exception
		{
		// Get size and rank in position communicator.
		p_size = p_comm.size();
		p_rank = p_comm.rank();

		// Create pseudorandom number generator.
		Random prng = Random.getInstance (seed);

		// Set up array slices.
		ranges = new Range (0, N-1) .subranges (p_size);
		myrange = ranges[p_rank];
		mylb = myrange.lb();
		myub = myrange.ub();
		mylength = myrange.length();

		// Initialize force, velocity, and position (x,y) arrays.
		f = new Vector2D [mylength];
		v = new Vector2D [mylength];
		for (int i = 0; i < mylength; ++ i)
			{
			f[i] = new Vector2D();
			v[i] = new Vector2D();
			}
		p = new Vector2D [N];
		for (int i = 0; i < N; ++ i)
			{
			p[i] =
				new Vector2D
					((0.9 * prng.nextDouble() + 0.05) * R,
					 (0.9 * prng.nextDouble() + 0.05) * R);
			}

		// Set up communication buffers.
		p_slices = Vector2D.doubleSliceBuffers (p, ranges);
		myp_slice = p_slices[p_rank];

		// Set up parallel team and parallel region.
		team = new ParallelTeam();
		region = new ParallelRegion()
			{
			public void run() throws Exception
				{
				// Calculate the net force on each antiproton.
				execute (mylb, myub, new IntegerForLoop()
					{
					public void run (int first, int last) throws Exception
						{
						double d;
						Vector2D temp = new Vector2D();

						// Accumulate forces between each pair of antiprotons,
						// but not between an antiproton and itself.
						for (int i = first; i <= last; ++ i)
							{
							Vector2D f_i = f[i-mylb];
							Vector2D p_i = p[i];
							for (int j = 0; j < i; ++ j)
								{
								Vector2D p_j = p[j];
								temp.assign (p_i);
								temp.sub (p_j);
								d = temp.mag();
								temp.mul (QP_QP / (d*d*d));
								f_i.add (temp);
								}
							for (int j = i+1; j < N; ++ j)
								{
								Vector2D p_j = p[j];
								temp.assign (p_i);
								temp.sub (p_j);
								d = temp.mag();
								temp.mul (QP_QP / (d*d*d));
								f_i.add (temp);
								}
							}

						for (int i = first; i <= last; ++ i)
							{
							Vector2D f_i = f[i-mylb];
							Vector2D p_i = p[i];

							// Accumulate force on antiproton from bottom side
							// of trap.
							double xp = p_i.x;
							double xpsqr = xp * xp;
							double xpmR = p_i.x - R;
							double xpmRsqr = xpmR * xpmR;
							double ypmyt = p_i.y;
							double ypmytsqr = ypmyt * ypmyt;
							f_i.x +=
								QP_QT *
									(1.0 / Math.sqrt (xpmRsqr + ypmytsqr) -
									 1.0 / Math.sqrt (xpsqr + ypmytsqr));
							f_i.y +=
								QP_QT / ypmyt *
									(-xpmR / Math.sqrt (xpmRsqr + ypmytsqr) +
									 xp / Math.sqrt (xpsqr + ypmytsqr));

							// Accumulate force on antiproton from top side of
							// trap.
							ypmyt = p_i.y - R;
							ypmytsqr = ypmyt * ypmyt;
							f_i.x +=
								QP_QT *
									(1.0 / Math.sqrt (xpmRsqr + ypmytsqr) -
									 1.0 / Math.sqrt (xpsqr + ypmytsqr));
							f_i.y +=
								QP_QT / ypmyt *
									(-xpmR / Math.sqrt (xpmRsqr + ypmytsqr) +
									 xp / Math.sqrt (xpsqr + ypmytsqr));

							// Accumulate force on antiproton from left side of
							// trap.
							double yp = p_i.y;
							double ypsqr = yp * yp;
							double ypmR = p_i.y - R;
							double ypmRsqr = ypmR * ypmR;
							double xpmxt = p_i.x;
							double xpmxtsqr = xpmxt * xpmxt;
							f_i.y +=
								QP_QT *
									(1.0 / Math.sqrt (ypmRsqr + xpmxtsqr) -
									 1.0 / Math.sqrt (ypsqr + xpmxtsqr));
							f_i.x +=
								QP_QT / xpmxt *
									(-ypmR / Math.sqrt (ypmRsqr + xpmxtsqr) +
									 yp / Math.sqrt (ypsqr + xpmxtsqr));

							// Accumulate force on antiproton from right side of
							// trap.
							xpmxt = p_i.x - R;
							xpmxtsqr = xpmxt * xpmxt;
							f_i.y +=
								QP_QT *
									(1.0 / Math.sqrt (ypmRsqr + xpmxtsqr) -
									 1.0 / Math.sqrt (ypsqr + xpmxtsqr));
							f_i.x +=
								QP_QT / xpmxt *
									(-ypmR / Math.sqrt (ypmRsqr + xpmxtsqr) +
									 yp / Math.sqrt (ypsqr + xpmxtsqr));
							}
						}
					});

				// The threads all wait at a barrier at the end of the above
				// parallel for loop before proceeding. This ensures that all
				// force calculations based on the antiprotons' positions for
				// the current time step have finished, before updating the
				// antiprotons' positions for the next time step.

				// Update each antiproton's position and velocity based on the
				// net force. No barrier wait needed at the end of this parallel
				// for loop, barrier wait at the end of the parallel region
				// suffices.
				execute (mylb, myub, new IntegerForLoop()
					{
					public void run (int first, int last) throws Exception
						{
						Vector2D temp = new Vector2D();

						for (int i = first; i <= last; ++ i)
							{
							Vector2D f_i = f[i-mylb];
							Vector2D v_i = v[i-mylb];
							Vector2D p_i = p[i];

							// Update antiproton's position and velocity.
							temp.assign (v_i);
							p_i.add (temp.mul (dt));
							temp.assign (f_i);
							p_i.add (temp.mul (onehalfdtsqr));
							temp.assign (f_i);
							v_i.add (temp.mul (dt));

							// Clear antiproton's force for the next
							// accumulation.
							f_i.clear();
							}
						}
					},
				BarrierAction.NO_WAIT);
				}
			};

		// Gather initial positions into process size-1.
		world.gather (size-1, myp_slice, null);

		// Do time steps.
		for (int frame = 1; frame <= frames; ++ frame)
			{
			for (int t = 0; t < steps; ++ t)
				{
				// Compute new antiproton positions in parallel threads.
				team.execute (region);

				// All-gather position array slices among processes 0 .. size-2.
				p_comm.allGather (myp_slice, p_slices);
				}

			// Gather positions after time steps into process size-1.
			world.gather (size-1, myp_slice, null);
			}
		}

	/**
	 * Visualization computations in process size-1.
	 */
	private static void computeVisualizations()
		throws Exception
		{
		// Set up position array.
		p = new Vector2D [N];
		for (int i = 0; i < N; ++ i)
			{
			p[i] = new Vector2D();
			}

		// Set up array slices.
		ranges = new Range [size];
		Range n_range = new Range (0, N-1);
		for (int i = 0; i <= size-2; ++ i)
			{
			ranges[i] = n_range.subrange (size-1, i);
			}
		ranges[size-1] = new Range(); // Zero-length range

		// Set up communication buffers.
		p_slices = Vector2D.doubleSliceBuffers (p, ranges);
		myp_slice = p_slices[size-1];

		// Set up for drawing.
		scale = ((double) W) / ((double) R);
		dot = new Ellipse2D.Double();
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

		// Render visualization frames.
		for (int frame = 0; frame <= frames; ++ frame)
			{
			// Gather antiproton positions from processes 0 .. size-2.
			world.gather (size-1, myp_slice, p_slices);

			// Fill in the background.
			g2d.setColor (BACKGROUND_COLOR);
			g2d.fillRect (0, 0, W, W);

			// Draw antiprotons.
			g2d.setColor (ANTIPROTON_COLOR);
			for (int i = 0; i < N; ++ i)
				{
				dot.setFrame
					(scale * p[i].x - DIAM_OVER_2,
					 scale * p[i].y - DIAM_OVER_2,
					 DIAM, DIAM);
				g2d.fill (dot);
				}

			// Generate file name.
			filename.setLength (0);
			filename.append (frame);
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
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<Kp> -Dpj.nt=<Kt> edu.rit.hyb.antimatter.AntiprotonHyb <seed> <N> <R> <frames> <steps> <dt> <W> <file>");
		System.err.println ("<Kp> = Number of parallel processes (>= 2)");
		System.err.println ("<Kt> = Number of parallel threads per process");
		System.err.println ("<seed> = Random seed");
		System.err.println ("<N> = Number of antiprotons");
		System.err.println ("<R> = Size of antiproton trap");
		System.err.println ("<frames> = Number of visualization frames");
		System.err.println ("<steps> = Number of time steps per frame");
		System.err.println ("<dt> = Size of time step");
		System.err.println ("<W> = Size of frame, <W>x<W> pixels");
		System.err.println ("<file> = Frame file names: \"<file>_0000.png\", etc.");
		System.exit (1);
		}

	}
