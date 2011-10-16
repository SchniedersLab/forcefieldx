//******************************************************************************
//
// File:    HotSpotClu.java
// Package: edu.rit.clu.heat
// Unit:    Class edu.rit.clu.heat.HotSpotClu
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

import edu.rit.io.Files;

import edu.rit.mp.DoubleBuf;

import edu.rit.mp.buf.DoubleItemBuf;

import edu.rit.pj.Comm;

import edu.rit.pj.reduction.DoubleOp;

import edu.rit.util.Arrays;
import edu.rit.util.Range;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;

/**
 * Class HotSpotClu is a cluster parallel program that calculates the
 * temperature distribution over a metal plate with hot spots.
 * <P>
 * Usage: java -Dpj.np=<I>K</I> edu.rit.clu.heat.HotSpotClu <I>imagefile</I>
 * <I>H</I> <I>W</I> <I>rl1</I> <I>cl1</I> <I>ru1</I> <I>cu1</I> <I>temp1</I> [
 * <I>rl2</I> <I>cl2</I> <I>ru2</I> <I>cu2</I> <I>temp2</I> . . . ]
 * <BR><I>K</I> = Number of parallel processes
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
 * red; intermediate temperatures are intermediate hues. Each process writes its
 * own PJG image file with a slice of the image. If <I>imagefile</I> is
 * specified as <TT>"out.pjg"</TT>, for example, then process 0 writes file
 * <TT>"out_0.pjg"</TT>, process 1 writes file <TT>"out_1.pjg"</TT>, and so on.
 * <P>
 * The computation is performed in parallel in multiple processors. The program
 * measures the computation's running time.
 *
 * @author  Alan Kaminsky
 * @version 11-Apr-2008
 */
public class HotSpotClu
	{

// Prevent construction.

	private HotSpotClu()
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

	private static final int FIRST  = 0;
	private static final int MIDDLE = 1;
	private static final int LAST   = 2;
	private static final int SINGLE = 3;

// Hidden variables.

	// World communicator.
	static Comm world;
	static int size;
	static int rank;
	static int position;
	static int predRank;
	static int succRank;

	// Command line arguments.
	static File imagefile;
	static int H;
	static int W;

	// Row slice index ranges.
	static Range[] slices;
	static Range mySlice;
	static int myLb;
	static int myUb;
	static int myLen;

	// Temperature grid.
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

	// Communication buffers.
	static DoubleItemBuf xibuf;
	static DoubleBuf hbuf_pred_red;
	static DoubleBuf hbuf_pred_black;
	static DoubleBuf hbuf_top_red;
	static DoubleBuf hbuf_top_black;
	static DoubleBuf hbuf_bottom_red;
	static DoubleBuf hbuf_bottom_black;
	static DoubleBuf hbuf_succ_red;
	static DoubleBuf hbuf_succ_black;

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
		if (size == 1) position = SINGLE;
		else if (rank == 0) position = FIRST;
		else if (rank < size-1) position = MIDDLE;
		else position = LAST;
		predRank = rank - 1;
		succRank = rank + 1;

		// Parse command line arguments.
		if (args.length < 8 || (args.length % 5) != 3) usage();
		imagefile = new File (args[0]);
		H = Integer.parseInt (args[1]);
		W = Integer.parseInt (args[2]);
		if (H < 1) usage();
		if (W < 1) usage();

		// Determine row slice index ranges.
		slices = new Range (1, H) .subranges (size);
		mySlice = slices[rank];
		myLb = mySlice.lb();
		myUb = mySlice.ub();
		myLen = mySlice.length();

		// Initialize temperature and hot spot meshes.
		h = new double [H+2] [];
		Arrays.allocate (h, new Range (myLb-1, myUb+1), W+2);
		hotspot = new boolean [H+2] [];
		Arrays.allocate (hotspot, new Range (myLb-1, myUb+1), W+2);

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
				if (h_r != null)
					{
					for (int c = cl; c <= cu; ++ c)
						{
						h_r[c] = temp;
						hotspot_r[c] = true;
						}
					}
				}
			}

		// Initialize communication buffers.
		xibuf = DoubleBuf.buffer();
		hbuf_pred_red = redBuffer (myLb-1);
		hbuf_pred_black = blackBuffer (myLb-1);
		hbuf_top_red = redBuffer (myLb);
		hbuf_top_black = blackBuffer (myLb);
		hbuf_bottom_red = redBuffer (myUb);
		hbuf_bottom_black = blackBuffer (myUb);
		hbuf_succ_red = redBuffer (myUb+1);
		hbuf_succ_black = blackBuffer (myUb+1);

		// Compute initial total absolute residual, then multiply by EPS.
		totalAbsXi = 0.0;
		double xi;
		for (int r = myLb; r <= myUb; ++ r)
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
		xibuf.item = totalAbsXi;
		world.allReduce (xibuf, DoubleOp.SUM);
		totalAbsXi = xibuf.item;
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
			for (int r = myLb; r <= myUb; ++ r)
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

			// Exchange boundary row red cells with neighboring processes.
			switch (position)
				{
				case FIRST:
					world.send (succRank, hbuf_bottom_red);
					world.receive (succRank, hbuf_succ_red);
					break;
				case MIDDLE:
					world.sendReceive
						(succRank, hbuf_bottom_red,
						 predRank, hbuf_pred_red);
					world.sendReceive
						(predRank, hbuf_top_red,
						 succRank, hbuf_succ_red);
					break;
				case LAST:
					world.receive (predRank, hbuf_pred_red);
					world.send (predRank, hbuf_top_red);
					break;
				}

			// Black half-sweep.
			for (int r = myLb; r <= myUb; ++ r)
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

			// Exchange boundary row black cells with neighboring processes.
			switch (position)
				{
				case FIRST:
					world.send (succRank, hbuf_bottom_black);
					world.receive (succRank, hbuf_succ_black);
					break;
				case MIDDLE:
					world.sendReceive
						(succRank, hbuf_bottom_black,
						 predRank, hbuf_pred_black);
					world.sendReceive
						(predRank, hbuf_top_black,
						 succRank, hbuf_succ_black);
					break;
				case LAST:
					world.receive (predRank, hbuf_pred_black);
					world.send (predRank, hbuf_top_black);
					break;
				}

			// Determine total absolute residual from all processes.
			xibuf.item = totalAbsXi;
			world.allReduce (xibuf, DoubleOp.SUM);
			totalAbsXi = xibuf.item;

			++ iterations;
			}
		while (iterations < MAXITER && totalAbsXi >= EPS_initialTotalAbsXi);

		// Check for convergence.
		if (iterations == MAXITER)
			{
			System.err.println ("HotSpotClu: Did not converge");
			System.exit (1);
			}

		long t3 = System.currentTimeMillis();

		// Generate image.
		int[][] matrix = new int [H+2] [];
		int rlb = rank == 0 ? myLb-1 : myLb;
		int rub = rank == size-1 ? myUb+1 : myUb;
		Arrays.allocate (matrix, new Range (rlb, rub), W+2);
		ColorImageRow matrix_r = new ColorImageRow (matrix[rlb]);
		for (int r = rlb; r <= rub; ++ r)
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
					(new FileOutputStream
						(Files.fileForRank (imagefile, rank))));
		writer.writeRowSlice (new Range (rlb, rub));
		writer.close();

		// Stop timing.
		long t4 = System.currentTimeMillis();
		System.out.println (iterations + " iterations " + rank);
		System.out.println ((t2-t1) + " msec pre " + rank);
		System.out.println ((t3-t2) + " msec calc " + rank);
		System.out.println ((t4-t3) + " msec post " + rank);
		System.out.println ((t4-t1) + " msec total " + rank);
		}

// Hidden operations.

	/**
	 * Returns a communication buffer for the red columns of the given row of
	 * the h matrix.
	 *
	 * @param  r  Row index.
	 *
	 * @return  Communication buffer.
	 */
	private static DoubleBuf redBuffer
		(int r)
		{
		return DoubleBuf.sliceBuffer (h[r], new Range (1 + (r&1), W, 2));
		}

	/**
	 * Returns a communication buffer for the black columns of the given row of
	 * the h matrix.
	 *
	 * @param  r  Row index.
	 *
	 * @return  Communication buffer.
	 */
	private static DoubleBuf blackBuffer
		(int r)
		{
		return DoubleBuf.sliceBuffer (h[r], new Range (2 - (r&1), W, 2));
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<K> edu.rit.clu.heat.HotSpotClu <imagefile> <H> <W> <rl1> <cl1> <ru1> <cu1> <temp1> [ <rl2> <cl2> <ru2> <cu2> <temp2> . . . ]");
		System.err.println ("<K> = Number of parallel processes");
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
