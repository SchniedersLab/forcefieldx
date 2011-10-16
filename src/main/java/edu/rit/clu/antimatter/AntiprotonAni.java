//******************************************************************************
//
// File:    AntiprotonAni.java
// Package: edu.rit.clu.antimatter
// Unit:    Class edu.rit.clu.antimatter.AntiprotonAni
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

package edu.rit.clu.antimatter;

import edu.rit.util.Random;

import edu.rit.vector.Vector2D;

import javax.swing.JFrame;

/**
 * Class AntiprotonAni is a sequential program that calculates the positions of
 * the antiprotons as a function of time and displays an animation on the
 * screen. Each antiproton experiences a net repulsive force from all the other
 * antiprotons. Each antiproton also experiences a force due to a magnetic field
 * perpendicular to the plane in which the antiprotons move.
 * <P>
 * Usage: java edu.rit.clu.antimatter.AntiprotonAni <I>seed</I> <I>R</I>
 * <I>dt</I> <I>steps</I> <I>N</I>
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
 * Displays a GUI showing the antiprotons' positions in the square from (0,0) to
 * (<I>R,R</I>). The total momentum of all the antiprotons is also shown.
 * <LI>
 * Performs <I>steps</I> time steps and updates the GUI. Each time step is
 * <I>dt</I>.
 * <LI>
 * Repeats until the GUI is closed or the program is externally terminated.
 * </OL>
 * <P>
 * The GUI window looks like this:
 * <P>
 * <CENTER>
 * <IMG SRC="doc-files/AntiprotonAni.png">
 * </CENTER>
 *
 * @author  Alan Kaminsky
 * @version 04-Feb-2008
 */
public class AntiprotonAni
	{

// Prevent construction.

	private AntiprotonAni()
		{
		}

// Hidden constants.

	// Charge on an antiproton.
	static final double QP = 3.0;

	// Magnetic field strength.
	static final double B = 3.0;

// Hidden variables.

	// Command line arguments.
	static long seed;
	static double R;
	static double dt;
	static int steps;
	static int N;

	// Acceleration, velocity, and position (x,y) arrays.
	static Vector2D[] a;
	static Vector2D[] v;
	static Vector2D[] p;

	// Total momentum.
	static Vector2D totalMV;

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
		seed = Long.parseLong (args[0]);
		R = Double.parseDouble (args[1]);
		dt = Double.parseDouble (args[2]);
		steps = Integer.parseInt (args[3]);
		N = Integer.parseInt (args[4]);

		double onehalfdtsqr = 0.5 * dt * dt;
		double QP_QP = QP * QP;
		double QP_B = QP * B;
		Vector2D temp = new Vector2D();

		// Create pseudorandom number generator.
		Random prng = Random.getInstance (seed);

		// Initialize acceleration, velocity, and position (x,y) arrays.
		a = new Vector2D [N];
		v = new Vector2D [N];
		p = new Vector2D [N];
		for (int i = 0; i < N; ++ i)
			{
			a[i] = new Vector2D();
			v[i] = new Vector2D();
			p[i] =
				new Vector2D
					(prng.nextDouble()*R/2+R/4, prng.nextDouble()*R/2+R/4);
			}
		totalMV = new Vector2D();

		// Create GUI.
		TrapFrame frame = new TrapFrame ("AntiprotonAni", p, R);
		TrapPanel panel = frame.getTrapPanel();
		frame.setDefaultCloseOperation (JFrame.EXIT_ON_CLOSE);
		frame.setVisible (true);

		// Do time steps forever.
		for (;;)
			{
			// Compute total momentum.
			totalMV.clear();
			for (int i = 0; i < N; ++ i)
				{
				totalMV.add (v[i]);
				}

			// Update display.
			frame.setTotalMomentum (totalMV.mag());
			panel.step();

			// Advance time by <steps> steps.
			for (int t = 0; t < steps; ++ t)
				{
				// Accumulate forces between each pair of antiprotons, but not
				// between an antiproton and itself.
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

				// Move all antiprotons.
				for (int i = 0; i < N; ++ i)
					{
					Vector2D a_i = a[i];
					Vector2D v_i = v[i];
					Vector2D p_i = p[i];

					// Accumulate force on antiproton from magnetic field.
					temp.assign (v_i) .mul (QP_B) .rotate270();
					a_i.add (temp);

					// Update antiproton's position and velocity.
					temp.assign (v_i);
					p_i.add (temp.mul (dt));
					temp.assign (a_i);
					p_i.add (temp.mul (onehalfdtsqr));
					temp.assign (a_i);
					v_i.add (temp.mul (dt));

					// Clear antiproton's acceleration for the next step.
					a_i.clear();
					}
				}
			}
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.clu.antimatter.AntiprotonAni <seed> <R> <dt> <steps> <N>");
		System.err.println ("<seed> = Random seed for initial antiproton positions");
		System.err.println ("<R> = Side of square for initial antiproton positions");
		System.err.println ("<dt> = Time step size");
		System.err.println ("<steps> = Number of time steps between snapshots");
		System.err.println ("<N> = Number of antiprotons");
		System.exit (1);
		}

	}
