//******************************************************************************
//
// File:    HotSpotSeq.java
// Package: edu.rit.clu.heat
// Unit:    Class edu.rit.clu.heat.HotSpotSeq
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

package edu.rit.clu.heat;

import edu.rit.color.HSB;

import edu.rit.image.ColorImageRow;
import edu.rit.image.PJGHueImage;
import edu.rit.image.PJGImage;

import edu.rit.pj.Comm;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;

/**
 * Class HotSpotSeq is a sequential program that calculates the temperature
 * distribution over a metal plate with hot spots.
 * <P>
 * Usage: java edu.rit.clu.heat.HotSpotSeq <I>imagefile</I> <I>H</I> <I>W</I>
 * <I>rl1</I> <I>cl1</I> <I>ru1</I> <I>cu1</I> <I>temp1</I>
 * [ <I>rl2</I> <I>cl2</I> <I>ru2</I> <I>cu2</I> <I>temp2</I> . . . ]
 * <BR><I>imagefile</I> = Output image file name
 * <BR><I>H</I> = Mesh height in pixels (<I>H</I> &gt;= 1)
 * <BR><I>W</I> = Mesh width in pixels (<I>W</I> &gt;= 1)
 * <BR><I>rl1</I> = First hot spot lower row (1 &lt;= <I>rl1</I> &lt;=
 * <I>W</I>)
 * <BR><I>cl1</I> = First hot spot lower column (1 &lt;= <I>cl1</I> &lt;=
 * <I>W</I>)
 * <BR><I>ru1</I> = First hot spot upper row (1 &lt;= <I>ru1</I> &lt;=
 * <I>W</I>)
 * <BR><I>cu1</I> = First hot spot upper column (1 &lt;= <I>cu1</I> &lt;=
 * <I>W</I>)
 * <BR><I>temp1</I> = First hot spot temperature (0.0 &lt;= <I>temp1</I> &lt;=
 * 100.0)
 * <P>
 * The program sets up a mesh of equally-spaced points with <I>H</I>+2 rows and
 * <I>W</I>+2 columns. The temperature of each boundary point [<I>r,c</I>],
 * where <I>r</I> = 0, <I>r</I> = <I>H</I>+1, <I>c</I> = 0, or <I>c</I> =
 * <I>W</I>+1, is fixed at 0.0. The temperatures at certain interior points,
 * known as "hot spots," are fixed at certain values greater than 0.0.
 * Specifically, the temperature of each point in the rectangle from
 * [<I>rl1,cl1</I>] to [<I>ru1,cu1</I>] inclusive is fixed at <I>temp1</I>, the
 * temperature of each point in the rectangle from [<I>rl2,cl2</I>] to
 * [<I>ru2,cu2</I>] inclusive is fixed at <I>temp2</I>, and so on. The program
 * calculates the temperature at every interior point (other than the hot spots)
 * using successive overrelaxation with Chebyshev acceleration and red-black
 * updating. The program outputs a PJG color image <I>H</I>+2 pixels high and
 * <I>W</I>+2 pixels wide. Each pixel's hue depends on the corresponding mesh
 * point's temperature. A temperature of 0.0 is blue; a temperature of 100.0 is
 * red; intermediate temperatures are intermediate hues. The PJG image file is
 * named <I>imagefile</I>.
 * <P>
 * The computation is performed sequentially in a single processor. The program
 * measures the computation's running time. This establishes a benchmark for
 * measuring the computation's running time on a parallel processor.
 *
 * @author  Alan Kaminsky
 * @version 11-Apr-2008
 */
public class HotSpotSeq
	{

// Prevent construction.

	private HotSpotSeq()
		{
		}

// Hidden constants.

	private static final double MIN_TEMP = 0.0;
	private static final double MAX_TEMP = 100.0;
	private static final double DELTA_TEMP = MAX_TEMP - MIN_TEMP;

	private static final double MIN_HUE = 4.0/6.0;
	private static final double MAX_HUE = 0.0;
	private static final double DELTA_HUE = MAX_HUE - MIN_HUE;

	private static final double EPS = 1.0e-3;

// Hidden variables.

	// Command line arguments.
	static File imagefile;
	static int H;
	static int W;

	// Temperature mesh.
	static double[][] h;

	// Mesh of hot spot locations.
	static boolean[][] hotspot;

	// Variables for total absolute residual.
	static double EPS_initialTotalAbsXi;
	static double totalAbsXi;

	// Other variables used in the successive overrrelaxation algorithm.
	static int MAXITER;
	static double rho_s_sqr;
	static double omega_over_4;
	static int iterations;

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
		if (args.length < 8 || (args.length % 5) != 3) usage();
		imagefile = new File (args[0]);
		H = Integer.parseInt (args[1]);
		W = Integer.parseInt (args[2]);
		if (H < 1) usage();
		if (W < 1) usage();

		// Initialize temperature and hot spot meshes.
		h = new double [H+2] [W+2];
		hotspot = new boolean [H+2] [W+2];

		// Record hot spot coordinates and temperatures.
		int n = (args.length - 3) / 5;
		for (int i = 0; i < n; ++ i)
			{
			int rl = Integer.parseInt (args[3+5*i]);
			int cl = Integer.parseInt (args[4+5*i]);
			int ru = Integer.parseInt (args[5+5*i]);
			int cu = Integer.parseInt (args[6+5*i]);
			double temp = Double.parseDouble (args[7+5*i]);
			if (1 > rl || rl > W) usage();
			if (1 > cl || cl > H) usage();
			if (1 > ru || ru > W) usage();
			if (1 > cu || cu > H) usage();
			if (MIN_TEMP > temp || temp > MAX_TEMP) usage();
			for (int r = rl; r <= ru; ++ r)
				{
				double[] h_r = h[r];
				boolean[] hotspot_r = hotspot[r];
				for (int c = cl; c <= cu; ++ c)
					{
					h_r[c] = temp;
					hotspot_r[c] = true;
					}
				}
			}

		// Compute initial total absolute residual, then multiply by EPS.
		totalAbsXi = 0.0;
		double xi;
		for (int r = 1; r <= H; ++ r)
			{
			double[] h_rm1 = h[r-1];
			double[] h_r   = h[r];
			double[] h_rp1 = h[r+1];
			boolean[] hotspot_r = hotspot[r];
			for (int c = 1; c <= W; ++ c)
				{
				xi =
					hotspot_r[c] ?
						0.0 :
						h_rm1[c]+h_rp1[c]+h_r[c-1]+h_r[c+1]-4.0*h_r[c];
				totalAbsXi += Math.abs (xi);
				}
			}
		EPS_initialTotalAbsXi = EPS * totalAbsXi;
		
		// Initialize other variables.
		MAXITER = 2 * (W + H);
		rho_s_sqr = 0.5 * (Math.cos (Math.PI/W) + Math.cos (Math.PI/H));
		rho_s_sqr = rho_s_sqr * rho_s_sqr;
		omega_over_4 = 0.25;
		iterations = 0;

		long t2 = System.currentTimeMillis();

		// Perform successive overrelaxation.
		do
			{
			totalAbsXi = 0.0;

			// Red half-sweep.
			for (int r = 1; r <= H; ++ r)
				{
				double[] h_rm1 = h[r-1];
				double[] h_r   = h[r];
				double[] h_rp1 = h[r+1];
				boolean[] hotspot_r = hotspot[r];
				for (int c = 1 + (r&1); c <= W; c += 2)
					{
					xi =
						hotspot_r[c] ?
							0.0 :
							h_rm1[c]+h_rp1[c]+h_r[c-1]+h_r[c+1]-4.0*h_r[c];
					totalAbsXi += Math.abs (xi);
					h_r[c] += omega_over_4 * xi;
					}
				}
			omega_over_4 = 0.25 /
				(1.0 - rho_s_sqr * (iterations == 0 ? 0.5 : omega_over_4));

			// Black half-sweep.
			for (int r = 1; r <= H; ++ r)
				{
				double[] h_rm1 = h[r-1];
				double[] h_r   = h[r];
				double[] h_rp1 = h[r+1];
				boolean[] hotspot_r = hotspot[r];
				for (int c = 2 - (r&1); c <= W; c += 2)
					{
					xi =
						hotspot_r[c] ?
							0.0 :
							h_rm1[c]+h_rp1[c]+h_r[c-1]+h_r[c+1]-4.0*h_r[c];
					totalAbsXi += Math.abs (xi);
					h_r[c] += omega_over_4 * xi;
					}
				}
			omega_over_4 = 0.25 / (1.0 - rho_s_sqr * omega_over_4);

			++ iterations;
			}
		while (iterations < MAXITER && totalAbsXi >= EPS_initialTotalAbsXi);

		// Check for convergence.
		if (iterations == MAXITER)
			{
			System.err.println ("HotSpotSeq: Did not converge");
			System.exit (1);
			}

		long t3 = System.currentTimeMillis();

		// Generate image.
		int[][] matrix = new int [H+2] [W+2];
		ColorImageRow matrix_r = new ColorImageRow (matrix[0]);
		for (int r = 0; r <= H+1; ++ r)
			{
			double[] h_r = h[r];
			matrix_r.setArray (matrix[r]);
			for (int c = 0; c <= W+1; ++ c)
				{
				matrix_r.setPixelHSB
					(/*c  */ c,
					 /*hue*/ (float)
						((h_r[c]-MIN_TEMP)/DELTA_TEMP*DELTA_HUE+MIN_HUE),
					 /*sat*/ 1.0f,
					 /*bri*/ 1.0f);
				}
			}
		PJGHueImage image = new PJGHueImage (H+2, W+2, matrix);
		PJGImage.Writer writer =
			image.prepareToWrite
				(new BufferedOutputStream
					(new FileOutputStream (imagefile)));
		writer.write();
		writer.close();

		// Stop timing.
		long t4 = System.currentTimeMillis();
		System.out.println (iterations + " iterations");
		System.out.println ((t2-t1) + " msec pre");
		System.out.println ((t3-t2) + " msec calc");
		System.out.println ((t4-t3) + " msec post");
		System.out.println ((t4-t1) + " msec total");
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.clu.heat.HotSpotSeq <imagefile> <H> <W> <rl1> <cl1> <ru1> <cu1> <temp1> [ <rl2> <cl2> <ru2> <cu2> <temp2> . . . ]");
		System.err.println ("<imagefile> = Output image file name");
		System.err.println ("<H> = Mesh height in pixels (<H> >= 1)");
		System.err.println ("<W> = Mesh width in pixels (<W> >= 1)");
		System.err.println ("<rl1> = First hot spot lower row (1 <= <rl1> <= <W>)");
		System.err.println ("<cl1> = First hot spot lower column (1 <= <cl1> <= <H>)");
		System.err.println ("<ru1> = First hot spot upper row (1 <= <ru1> <= <W>)");
		System.err.println ("<cu1> = First hot spot upper column (1 <= <cu1> <= <H>)");
		System.err.println ("<temp1> = First hot spot temperature (0.0 <= <temp1> <= 100.0)");
		System.exit (1);
		}

	}
