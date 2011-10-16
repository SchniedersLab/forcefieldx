//******************************************************************************
//
// File:    PixelAnalysis.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.PixelAnalysis
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

import edu.rit.mri.SpinSignal;
import edu.rit.mri.SpinSignalDifference;

import edu.rit.numeric.NonLinearLeastSquares;
import edu.rit.numeric.NonNegativeLeastSquares;
import edu.rit.numeric.Series;
import edu.rit.numeric.TooManyIterationsException;

import java.util.ArrayList;
import java.util.List;

/**
 * Class PixelAnalysis provides a routine for doing a spin relaxometry analysis
 * on one pixel of a magnetic resonance image.
 * <P>
 * The input to the analysis is a measured spin signal expressed as two
 * {@linkplain edu.rit.numeric.Series Series} objects, a time series
 * <I>t</I><SUB><I>i</I></SUB> and a spin signal series
 * <I>S</I>(<I>t</I><SUB><I>i</I></SUB>). Another input is a series of fixed
 * spin-lattice relaxation rates <I>R</I>1<SUB><I>j</I></SUB>. These rates are
 * chosen to cover the range of likely rates for the magnetic resonance image
 * being analyzed.
 * <P>
 * The routine first does a <I>nonnegative, linear</I> least squares fit of the
 * input data to a model consisting of a group of tissues with the input fixed
 * spin-lattice relaxation rates <I>R</I>1<SUB><I>j</I></SUB>. Peaks in the
 * linear least squares fit determine the number of tissues and the approximate
 * spin density and spin-lattice relaxation rate for each tissue.
 * <P>
 * The routine then does a <I>nonlinear</I> least squares fit of the input data
 * to a model consisting of the number of tissues determined in the previous
 * step. The nonlinear least squares fit "polishes up" the spin densities and
 * spin-lattice relaxation rates determined in the previous step, which are only
 * approximate.
 * <P>
 * The routine checks the nonlinear least squares fit for plausibility. To be
 * plausible:
 * <UL>
 * <LI>
 * All spin densities and spin-lattice relaxation rates must be positive. (If
 * this is not the case, it's likely the nonlinear least squares fit is trying
 * to fit the data to too many parameters, resulting in nonsensical parameter
 * values.)
 * <P><LI>
 * All the spin-lattice relaxation rates must be sufficiently far apart.
 * Specifically, the relative difference between any two rates must be greater
 * than 0.001. (If the rates are closer together than that, it's likely they
 * represent the same tissue.)
 * <P><LI>
 * The sum of the spin densities must agree with the asymptotic spin signal for
 * large values of <I>t</I>. Specifically, the sum must be within 20% of the
 * average of the last seven spin signal values. (If this is not the case, it's
 * likely that two spurious "tissues," one with a very large relaxation rate,
 * one with a very small relaxation rate, are canceling each other out, and
 * there really should be only one tissue.)
 * </UL>
 * If the nonlinear least squares fit is not plausible, the routine decides it
 * is trying to fit too many tissues. The routine eliminates the tissue with the
 * smallest spin density and repeats the nonlinear least squares fit. This
 * continues until the fit is plausible or until all the tissues have been
 * eliminated, in which case the routine reports that it could not find a
 * solution.
 * <P>
 * The output of the analysis is a list of the tissues' computed spin densities
 * and a list of the tissues' computed spin-lattice relaxation rates.
 *
 * @author  Alan Kaminsky
 * @version 16-Jun-2008
 */
public class PixelAnalysis
	{

// Prevent construction.

	private PixelAnalysis()
		{
		}

// Exported operations.

	/**
	 * Do a spin relaxometry analysis.
	 *
	 * @param  t_series
	 *     Series of measured time values, of length <I>M</I> (input).
	 * @param  S_series
	 *     Series of measured spin signal values, of length <I>M</I> (input).
	 * @param  R1_series
	 *     Series of fixed spin-lattice relaxation rates for the linear part of
	 *     the analysis, of length <I>N</I> (input).
	 * @param  A
	 *     Design matrix for the linear part of the analysis (input). This must
	 *     be an <I>M</I>&times;<I>N</I>-element matrix such that
	 *     <I>A</I><SUB><I>i,j</I></SUB> = 1 &minus; 2
	 *     exp(&minus;<I>R</I>1<SUB><I>j</I></SUB> <I>t</I><SUB><I>i</I></SUB>).
	 *     (The design matrix is supplied as an argument because the same design
	 *     matrix is typically used for every pixel in an image, and calculating
	 *     the design matrix just once outside this routine saves time.)
	 * @param  rho_list
	 *     List in which to store the computed spin densities (output). The size
	 *     of the list is the number of tissues. If the routine could not find a
	 *     solution, the size of the list is 0.
	 * @param  R1_list
	 *     List in which to store the computed spin-lattice relaxation rates
	 *     (output). The size of the list is the number of tissues. If the
	 *     routine could not find a solution, the size of the list is 0.
	 */
	public static void analyze
		(Series t_series,
		 Series S_series,
		 Series R1_series,
		 double[][] A,
		 List<Double> rho_list,
		 List<Double> R1_list)
		{
		int M = t_series.length();
		int N = R1_series.length();

		// Do a spin relaxometry analysis using nonnegative linear least
		// squares.

		// Create nonnegative linear least squares solver.
		NonNegativeLeastSquares linsolver = new NonNegativeLeastSquares (M, N);

		// Find the solution.
		for (int i = 0; i < M; ++ i)
			{
			System.arraycopy (A[i], 0, linsolver.a[i], 0, N);
			linsolver.b[i] = S_series.x(i);
			}
		linsolver.solve();
		double[] rho_series = linsolver.x;

		// Find peaks in the solution. A peak occurs at index i if
		// rho[i] > rho[i-1] and rho[i] > rho[i+1].
		ArrayList<Double> approx_rho_list = new ArrayList<Double>();
		ArrayList<Double> approx_R1_list = new ArrayList<Double>();
		for (int j = 0; j < N; ++ j)
			{
			if (rho_series[j] > (j == 0 ? 0.0 : rho_series[j-1]) &&
					rho_series[j] > (j == N-1 ? 0.0 : rho_series[j+1]))
				{
				approx_rho_list.add (rho_series[j]);
				approx_R1_list.add (R1_series.x(j));
				}
			}

		// Do a spin relaxometry analysis using nonlinear least squares. Peaks
		// in the linear analysis give the initial vector of densities and
		// rates.

		// Repeat until the solution is plausible.
		boolean plausible = false;
		int L = approx_rho_list.size();
		rho_list.clear();
		R1_list.clear();
		while (L > 0 && ! plausible)
			{
			// Create spin signal difference function. L = number of tissues.
			SpinSignalDifference fcn =
				new SpinSignalDifference (t_series, S_series, L);

			// Create nonlinear least squares solver.
			NonLinearLeastSquares nonlinsolver =
				new NonLinearLeastSquares (fcn);

			// Find the solution.
			for (int i = 0; i < L; ++ i)
				{
				nonlinsolver.x[(i<<1)] = approx_rho_list.get(i);
				nonlinsolver.x[(i<<1)+1] = approx_R1_list.get(i);
				}
			try
				{
				nonlinsolver.solve();
				for (int i = 0; i < L; ++ i)
					{
					rho_list.add (nonlinsolver.x[(i<<1)]);
					R1_list.add (nonlinsolver.x[(i<<1)+1]);
					}

				// Decide if solution is plausible.
				plausible = checkPlausibility (S_series, rho_list, R1_list);
				}

			// Couldn't find a solution.
			catch (TooManyIterationsException exc)
				{
				plausible = false;
				}

			// If solution is not plausible, eliminate tissue with smallest
			// density and try again.
			if (! plausible)
				{
				double minrho = Double.MAX_VALUE;
				int mini = 0;
				for (int i = 0; i < L; ++ i)
					{
					if (approx_rho_list.get(i) < minrho)
						{
						minrho = approx_rho_list.get(i);
						mini = i;
						}
					}
				approx_rho_list.remove (mini);
				approx_R1_list.remove (mini);
				L = approx_rho_list.size();
				rho_list.clear();
				R1_list.clear();
				}
			}
		}

// Hidden operations.

	/**
	 * Decide if the given solution is plausible.
	 */
	private static boolean checkPlausibility
		(Series S_series,
		 List<Double> rho_list,
		 List<Double> R1_list)
		{
		int M = S_series.length();
		int L = rho_list.size();

		// If any density or rate is negative, solution is not plausible.
		for (int i = 0; i < L; ++ i)
			{
			if (rho_list.get(i) < 0.0)
				{
				return false;
				}
			if (R1_list.get(i) < 0.0)
				{
				return false;
				}
			}

		// If relative difference between any two rates is too small, solution
		// is not plausible.
		for (int i = 0; i < L-1; ++ i)
			{
			double R_i = R1_list.get(i);
			for (int j = i+1; j < L; ++ j)
				{
				double R_j = R1_list.get(j);
				double reldiff = 2.0*Math.abs(R_i-R_j)/Math.abs(R_i+R_j);
				if (reldiff <= 0.001)
					{
					return false;
					}
				}
			}

		// If sum of densities is too far from asymptotic measurement for large
		// t, solution is not plausible.
		double sumrho = 0.0;
		for (int i = 0; i < L; ++ i)
			{
			sumrho += rho_list.get(i);
			}
		double S_last = 0.0;
		int n = 0;
		for (int i = M-1; i >=0 && n < 7; -- i)
			{
			S_last += S_series.x(i);
			++ n;
			}
		S_last /= n;
		double reldiff = Math.abs(sumrho-S_last)/Math.abs(S_last);
		if (reldiff >= 0.2)
			{
			return false;
			}

		// Solution is plausible.
		return true;
		}

	}
