//******************************************************************************
//
// File:    AntiprotonClu.java
// Package: edu.rit.clu.antimatter
// Unit:    Class edu.rit.clu.antimatter.AntiprotonClu
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

import edu.rit.io.Files;

import edu.rit.mp.DoubleBuf;

import edu.rit.pj.Comm;

import edu.rit.util.Random;
import edu.rit.util.Range;

import edu.rit.vector.Vector2D;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;

/**
 * Class AntiprotonClu is a cluster parallel program that calculates the
 * positions of the antiprotons as a function of time. Each antiproton
 * experiences a net repulsive force from all the other antiprotons. Each
 * antiproton also experiences a force due to a magnetic field perpendicular to
 * the plane in which the antiprotons move.
 * <P>
 * Usage: java -Dpj.np=<I>K</I> edu.rit.clu.antimatter.AntiprotonClu <I>seed</I>
 * <I>R</I> <I>dt</I> <I>steps</I> <I>snaps</I> <I>N</I> <I>outfile</I>
 * <P>
 * The program runs in <I>K</I> parallel processes. Each process writes its own
 * output file. If <I>outfile</I> is specified as <TT>"out.dat"</TT>, for
 * example, then process 0 writes file <TT>"out_0.dat"</TT>, process 1 writes
 * file <TT>"out_1.dat"</TT>, and so on.
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
 * The computation is performed in parallel in multiple processors. The program
 * measures the computation's running time.
 * <P>
 * Each process holds the entire position array, one slice of the velocity
 * array, and one slice of the acceleration array. At each time step, each
 * process calculates its own slice of the acceleration array using the entire
 * position array. Each process uses its slice of the acceleration array to
 * update its slices of the velocity and position arrays. Then the processes do
 * an all-gather of the position array slices, so every processor has the
 * complete new position array for the next time step. Each process writes
 * snapshots of its own slice of the position array into the process's own
 * output file.
 *
 * @author  Alan Kaminsky
 * @version 09-Feb-2008
 */
public class AntiprotonClu
	{

// Prevent construction.

	private AntiprotonClu()
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

	// World communicator.
	static Comm world;
	static int size;
	static int rank;

	// Command line arguments.
	static long seed;
	static double R;
	static double dt;
	static int steps;
	static int snaps;
	static int N;
	static File outfile;

	static double one_half_dt_sqr;

	// Antiproton slices.
	static Range[] slices;
	static Range mySlice;
	static int myLb;
	static int myLen;

	// Acceleration, velocity, and position vector arrays.
	static Vector2D[] a;
	static Vector2D[] v;
	static Vector2D[] p;

	// Total momentum.
	static Vector2D totalMV = new Vector2D();

	// Position array communication buffers.
	static DoubleBuf[] buffers;
	static DoubleBuf myBuffer;

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
		// Start timing.
		long t1 = System.currentTimeMillis();

		// Initialize world communicator.
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();

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

		// Set up antiproton slices.
		slices = new Range (0, N-1) .subranges (size);
		mySlice = slices[rank];
		myLb = mySlice.lb();
		myLen = mySlice.length();

		// Create pseudorandom number generator.
		Random prng = Random.getInstance (seed);

		// Initialize position vector array with all antiprotons.
		p = new Vector2D [N];
		for (int i = 0; i < N; ++ i)
			{
			p[i] = new Vector2D
				(prng.nextDouble()*R/2+R/4, prng.nextDouble()*R/2+R/4);
			}

		// Initialize acceleration and velocity vector arrays with a slice of
		// antiprotons.
		a = new Vector2D [myLen];
		v = new Vector2D [myLen];
		for (int i = 0; i < myLen; ++ i)
			{
			a[i] = new Vector2D();
			v[i] = new Vector2D();
			}

		// Set up position array communication buffers.
		buffers = Vector2D.doubleSliceBuffers (p, slices);
		myBuffer = buffers[rank];

		// Set up output file and write initial snapshot.
		AntiprotonFile out =
			new AntiprotonFile (seed, R, dt, steps, snaps+1, N, myLb, myLen);
		AntiprotonFile.Writer writer =
			out.prepareToWrite
				(new BufferedOutputStream
					(new FileOutputStream
						(Files.fileForRank (outfile, rank))));
		writer.writeSnapshot (p, myLb, totalMV);

		long t2 = System.currentTimeMillis();

		// Do <snaps> snapshots.
		for (int s = 0; s < snaps; ++ s)
			{
			// Advance time by <steps> steps.
			for (int t = 0; t < steps; ++ t)
				{
				computeAcceleration();
				step();

				// All-gather the new antiproton positions.
				world.allGather (myBuffer, buffers);
				}

			// Compute total momentum.
			computeTotalMomentum();

			// Write snapshot.
			writer.writeSnapshot (p, myLb, totalMV);
			}

		// Close output file.
		writer.close();

		// Stop timing.
		long t3 = System.currentTimeMillis();
		System.out.println ((t2-t1) + " msec pre " + rank);
		System.out.println ((t3-t2) + " msec calc " + rank);
		System.out.println ((t3-t1) + " msec total " + rank);
		}

// Hidden operations.

	/**
	 * Compute this process's slice of the antiproton accelerations due to the
	 * repulsive forces from all the antiprotons.
	 */
	private static void computeAcceleration()
		{
		// Accumulate forces between each pair of antiprotons, but not between
		// an antiproton and itself.
		for (int i = 0; i < myLen; ++ i)
			{
			Vector2D a_i = a[i];
			int index = i + myLb;
			Vector2D p_i = p[index];
			for (int j = 0; j < index; ++ j)
				{
				temp.assign (p_i);
				temp.sub (p[j]);
				double dsqr = temp.sqrMag();
				temp.mul (QP_QP / (dsqr * Math.sqrt(dsqr)));
				a_i.add (temp);
				}
			for (int j = index+1; j < N; ++ j)
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
		// Move all antiprotons in this slice.
		for (int i = 0; i < myLen; ++ i)
			{
			Vector2D a_i = a[i];
			Vector2D v_i = v[i];
			Vector2D p_i = p[i+myLb];

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
	 * Compute the total momentum for this process's slice of the antiprotons.
	 * The answer is stored in <TT>totalMV</TT>.
	 */
	private static void computeTotalMomentum()
		{
		totalMV.clear();
		for (int i = 0; i < myLen; ++ i)
			{
			totalMV.add (v[i]);
			}
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<K> edu.rit.clu.antimatter.AntiprotonClu <seed> <R> <dt> <steps> <snaps> <N> <outfile>");
		System.err.println ("<K> = Number of parallel processes");
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
