//******************************************************************************
//
// File:    AntiprotonSeq.java
// Package: edu.rit.clu.antimatter
// Unit:    Class edu.rit.clu.antimatter.AntiprotonSeq
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

package edu.rit.clu.antimatter;

import edu.rit.pj.Comm;

import edu.rit.util.Random;
import edu.rit.util.Range;

import edu.rit.vector.Vector2D;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;

/**
 * Class AntiprotonSeq is a sequential program that calculates the positions of
 * the antiprotons as a function of time. Each antiproton experiences a net
 * repulsive force from all the other antiprotons. Each antiproton also
 * experiences a force due to a magnetic field perpendicular to the plane in
 * which the antiprotons move.
 * <P>
 * Usage: java edu.rit.clu.antimatter.AntiprotonSeq <I>seed</I> <I>R</I>
 * <I>dt</I> <I>steps</I> <I>snaps</I> <I>N</I> <I>outfile</I>
 * <P>
 * The program:
 * <OL TYPE=1>
 * <LI>
 * Initializes a pseudorandom number generator with <I>seed</I>.
 * <LI>
 * Generates <I>N</I> antiprotons positioned at random in the square from
 * (0.25<I>R</I>,0.25<I>R</I>) to (0.75<I>R</I>,0.75<I>R</I>).
 * <LI>
 * Sets each antiproton's initial velocity to 0.
 * <LI>
 * Stores a snapshot of the antiprotons' initial positions in the
 * <I>outfile</I>.
 * <LI>
 * Performs <I>steps</I> time steps and stores another snapshot of the
 * antiprotons' positions in the <I>outfile</I>. Each time step is <I>dt</I>.
 * <LI>
 * Repeats Step 5 <I>snaps</I> times. The number of snapshots stored in the
 * <I>outfile</I> is <I>snaps</I>+1.
 * </OL>
 * <P>
 * The computation is performed sequentially in a single processor. The program
 * measures the computation's running time. This establishes a benchmark for
 * measuring the computation's running time on a parallel processor.
 *
 * @author  Alan Kaminsky
 * @version 09-Feb-2008
 */
public class AntiprotonSeq
	{

// Prevent construction.

	private AntiprotonSeq()
		{
		}

// Hidden constants.

	// Charge on an antiproton.
	static final double QP = 3.0;

	// Magnetic field strength.
	static final double B = 3.0;

	static final double QP_QP = QP * QP;
	static final double QP_B = QP * B;

// Hidden variables.

	// Command line arguments.
	static long seed;
	static double R;
	static double dt;
	static int steps;
	static int snaps;
	static int N;
	static File outfile;

	static double one_half_dt_sqr;

	// Acceleration, velocity, and position vector arrays.
	static Vector2D[] a;
	static Vector2D[] v;
	static Vector2D[] p;

	// Total momentum.
	static Vector2D totalMV = new Vector2D();

	// Temporary storage.
	static Vector2D temp = new Vector2D();

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
		if (args.length != 7) usage();
		seed = Long.parseLong (args[0]);
		R = Double.parseDouble (args[1]);
		dt = Double.parseDouble (args[2]);
		steps = Integer.parseInt (args[3]);
		snaps = Integer.parseInt (args[4]);
		N = Integer.parseInt (args[5]);
		outfile = new File (args[6]);

		one_half_dt_sqr = 0.5 * dt * dt;

		// Create pseudorandom number generator.
		Random prng = Random.getInstance (seed);

		// Initialize acceleration, velocity, and position vector arrays.
		a = new Vector2D [N];
		v = new Vector2D [N];
		p = new Vector2D [N];
		for (int i = 0; i < N; ++ i)
			{
			a[i] = new Vector2D();
			v[i] = new Vector2D();
			p[i] = new Vector2D
				(prng.nextDouble()*R/2+R/4, prng.nextDouble()*R/2+R/4);
			}

		// Set up output file and write initial snapshot.
		AntiprotonFile out =
			new AntiprotonFile (seed, R, dt, steps, snaps+1, N, 0, N);
		AntiprotonFile.Writer writer =
			out.prepareToWrite
				(new BufferedOutputStream
					(new FileOutputStream (outfile)));
		writer.writeSnapshot (p, 0, totalMV);

		long t2 = System.currentTimeMillis();

		// Do <snaps> snapshots.
		for (int s = 0; s < snaps; ++ s)
			{
			// Advance time by <steps> steps.
			for (int t = 0; t < steps; ++ t)
				{
				computeAcceleration();
				step();
				}

			// Compute total momentum.
			computeTotalMomentum();

			// Write snapshot.
			writer.writeSnapshot (p, 0, totalMV);
			}

		// Close output file.
		writer.close();

		// Stop timing.
		long t3 = System.currentTimeMillis();
		System.out.println ((t2-t1) + " msec pre");
		System.out.println ((t3-t2) + " msec calc");
		System.out.println ((t3-t1) + " msec total");
		}

// Hidden operations.

	/**
	 * Compute the antiproton accelerations due to the repulsive forces from all
	 * the antiprotons.
	 */
	private static void computeAcceleration()
		{
		// Accumulate forces between each pair of antiprotons, but not between
		// an antiproton and itself.
		for (int i = 0; i < N; ++ i)
			{
			Vector2D a_i = a[i];
			Vector2D p_i = p[i];
			for (int j = 0; j < i; ++ j)
				{
				temp.assign (p_i);
				temp.sub (p[j]);
				double dsqr = temp.sqrMag();
				temp.mul (QP_QP / (dsqr * Math.sqrt(dsqr)));
				a_i.add (temp);
				}
			for (int j = i+1; j < N; ++ j)
				{
				temp.assign (p_i);
				temp.sub (p[j]);
				double dsqr = temp.sqrMag();
				temp.mul (QP_QP / (dsqr * Math.sqrt(dsqr)));
				a_i.add (temp);
				}
			}
		}

	/**
	 * Take one time step.
	 */
	private static void step()
		{
		// Move all antiprotons.
		for (int i = 0; i < N; ++ i)
			{
			Vector2D a_i = a[i];
			Vector2D v_i = v[i];
			Vector2D p_i = p[i];

			// Accumulate acceleration on antiproton from magnetic field.
			temp.assign (v_i) .mul (QP_B) .rotate270();
			a_i.add (temp);

			// Update antiproton's position and velocity.
			temp.assign (v_i);
			p_i.add (temp.mul (dt));
			temp.assign (a_i);
			p_i.add (temp.mul (one_half_dt_sqr));
			temp.assign (a_i);
			v_i.add (temp.mul (dt));

			// Clear antiproton's acceleration for the next step.
			a_i.clear();
			}
		}

	/**
	 * Compute the total momentum for all the antiprotons. The answer is stored
	 * in <TT>totalMV</TT>.
	 */
	private static void computeTotalMomentum()
		{
		totalMV.clear();
		for (int i = 0; i < N; ++ i)
			{
			totalMV.add (v[i]);
			}
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.clu.antimatter.AntiprotonSeq <seed> <R> <dt> <steps> <snaps> <N> <outfile>");
		System.err.println ("<seed> = Random seed for initial antiproton positions");
		System.err.println ("<R> = Side of square for initial antiproton positions");
		System.err.println ("<dt> = Time step size");
		System.err.println ("<steps> = Number of time steps between snapshots");
		System.err.println ("<snaps> = Number of snapshots");
		System.err.println ("<N> = Number of antiprotons");
		System.err.println ("<outfile> = Output file name");
		System.exit (1);
		}

	}
