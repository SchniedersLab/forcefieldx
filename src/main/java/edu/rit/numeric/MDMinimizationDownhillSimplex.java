//******************************************************************************
//
// File:    MDMinimizationDownhillSimplex.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.MDMinimizationDownhillSimplex
//
// This Java source file is copyright (C) 2010 by Alan Kaminsky. All rights
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

package edu.rit.numeric;

import static java.lang.Math.*;

import java.util.Arrays;

/**
 * Class MDMinimizationDownhillSimplex finds a minimum of a multidimensional
 * function using the downhill simplex method of Nelder and Mead. The function
 * has <I>N</I> inputs and is represented by an object that implements interface
 * {@linkplain MDFunction}. The <TT>minimize()</TT> method finds a (local)
 * minimum in the function. The input to the <TT>minimize()</TT> method is an
 * <I>N</I>-dimensional <I>simplex,</I> namely a group of <I>N</I>+1 points in
 * <I>N</I> dimensions, each point different from all the others. Typically, one
 * point of the simplex is an initial guess for the solution, and the other
 * points of the simplex are perturbations of the initial guess. The
 * <TT>minimize()</TT> method evaluates the function at each point of the
 * simplex and moves the simplex through <I>N</I>-dimensional space to minimize
 * the smallest function value of any point in the simplex, stopping when the
 * relative difference between the smallest and largest function values of any
 * points in the simplex falls below a tolerance. The output from the
 * <TT>minimize()</TT> method is the final simplex, with the
 * smallest-function-value point at position 0. The inputs to and outputs from
 * the <TT>minimize()</TT> method are stored in the fields of an instance of
 * class MDMinimizationDownhillSimplex.
 *
 * @author  Alan Kaminsky
 * @version 05-Oct-2010
 */
public class MDMinimizationDownhillSimplex
	{

// Exported data members.

	/**
	 * The multidimensional function to be minimized.
	 */
	public final MDFunction fcn;

	/**
	 * The number of function arguments.
	 */
	public final int N;

	/**
	 * The simplex. This is an <I>N</I>+1-element array of <I>N</I>-element
	 * points. On input to the <TT>minimize()</TT> method, <TT>x</TT> contains a
	 * simplex with an initial estimate of the solution. On output from the
	 * <TT>minimize()</TT> method, <TT>x</TT> contains the final simplex, and
	 * <TT>x[0]</TT> contains the solution, namely the simplex point with the
	 * smallest function value.
	 */
	public final double[][] x;

	/**
	 * Function values of the simplex points. An <I>N</I>+1-element array. On
	 * output from the <TT>minimize()</TT> method, <TT>f[i]</TT> contains the
	 * function value for simplex point <TT>x[i]</TT>, 0 &le; <TT>i</TT> &le;
	 * <I>N</I>&minus;1. In particular, <TT>f[0]</TT> contains the minimized
	 * function value for the solution.
	 */
	public final double[] f;

	/**
	 * Tolerance. An input to the <TT>minimize()</TT> method. Must be &gt; 0.
	 * Termination occurs when the relative difference between the smallest and
	 * largest function values of points in the simplex is less than the
	 * tolerance. The default tolerance is 1&times;10<SUP>&minus;6</SUP>.
	 */
	public double tol = 1.0e-6;

	/**
	 * Debug flag. An input to the <TT>minimize()</TT> method. If true, the
	 * <TT>subclassDebug()</TT> method is called at the beginning of every
	 * iteration. The default setting is false.
	 */
	public boolean debug;

// Hidden data members.

	// Indexes of simplex points with largest, second-largest, and smallest
	// function values.
	private int i_max;
	private int i_2ndmax;
	private int i_min;

	// Number of function evaluations.
	private int evalCount;

	// Sum of the simplex points.
	private double[] x_sum;

	// Trial simplex point, and its function value.
	private double[] x_trial;
	private double f_trial;

	// Maximum number of function evaluations allowed.
	private static final int MAXEVAL = 5000;

// Exported constructors.

	/**
	 * Construct a new multidimensional minimization object for the given
	 * function. Field <TT>fcn</TT> is set to <TT>theFunction</TT>. Field
	 * <TT>N</TT> is set by calling the function's <TT>argumentLength()</TT>
	 * method. Field <TT>x</TT>, the simplex, is allocated with the proper size;
	 * all simplex points are initially 0. Field <TT>f</TT>, the function
	 * values, is allocated with the proper size.
	 *
	 * @param  theFunction  Multidimensional function to be minimized.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFunction</TT> is null.
	 */
	public MDMinimizationDownhillSimplex
		(MDFunction theFunction)
		{
		fcn = theFunction;
		N = theFunction.argumentLength();
		x = new double [N + 1] [N];
		f = new double [N + 1];
		x_sum = new double [N];
		x_trial = new double [N];
		}

// Exported operations.

	/**
	 * Set the simplex to the given point plus the same perturbation along each
	 * dimension. The field <TT>x[0]</TT> is set to <TT>x</TT>; the field
	 * <TT>x[1]</TT> is set to <TT>x</TT> with <TT>delta</TT> added to element
	 * 0; the field <TT>x[2]</TT> is set to <TT>x</TT> with <TT>delta</TT> added
	 * to element 1; and so on.
	 *
	 * @param  x      Simplex point. Must be an <I>N</I>-element array.
	 * @param  delta  Perturbation along every dimension. Must be nonzero.
	 */
	public void setSimplex
		(double[] x,
		 double delta)
		{
		System.arraycopy (x, 0, this.x[0], 0, N);
		for (int i = 0; i < N; ++ i) this.x[i+1][i] += delta;
		}

	/**
	 * Set the simplex to the given point plus a different perturbation along
	 * each dimension. The field <TT>x[0]</TT> is set to <TT>x</TT>; the field
	 * <TT>x[1]</TT> is set to <TT>x</TT> with <TT>delta[0]</TT> added to
	 * element 0; the field <TT>x[2]</TT> is set to <TT>x</TT> with
	 * <TT>delta[1]</TT> added to element 1; and so on.
	 *
	 * @param  x      Simplex point. Must be an <I>N</I>-element array.
	 * @param  delta  Perturbations along each dimension. All must be nonzero.
	 */
	public void setSimplex
		(double[] x,
		 double[] delta)
		{
		System.arraycopy (x, 0, this.x[0], 0, N);
		for (int i = 0; i < N; ++ i) this.x[i+1][i] += delta[i];
		}

	/**
	 * Minimize the multidimensional function. On input, the field <TT>x</TT>
	 * must be filled in with an initial simplex, and the field <TT>tol</TT>
	 * must be set to the desired tolerance. On output, the field <TT>x</TT>
	 * contains the solution. For further information, see the documentation for
	 * each field.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>tol</TT> &le; 0.
	 * @exception  TooManyIterationsException
	 *     (unchecked exception) Thrown if too many function evaluations (5,000)
	 *     occurred without finding a minimum.
	 */
	public void minimize()
		{
		// Verify preconditions.
		if (tol <= 0.0)
			{
			throw new IllegalArgumentException
				("MDMinimizationDownhillSimplex.minimize(): tol = "+tol+
				 " illegal");
			}

		// Compute function values for initial simplex.
		for (int i = 0; i <= N; ++ i)
			{
			f[i] = fcn.f (x[i]);
			}
		evalCount = N + 1;

		// Compute sum of initial simplex points.
		compute_x_sum();

		// Iterate until solution is found or too many function evaluations.
		for (int iter = 1; ; ++ iter)
			{
			if (debug) subclassDebug (iter, evalCount);

			// Find largest, second-largest, and smallest simplex points.
			if (f[0] > f[1])
				{
				i_max = 0;
				i_2ndmax = i_min = 1;
				}
			else
				{
				i_max = 1;
				i_2ndmax = i_min = 0;
				}
			for (int i = 2; i <= N; ++ i)
				{
				if (f[i] > f[i_max])
					{
					i_2ndmax = i_max;
					i_max = i;
					}
				else if (f[i] > f[i_2ndmax])
					{
					i_2ndmax = i;
					}
				if (f[i] < f[i_min])
					{
					i_min = i;
					}
				}

			// Compute relative difference between largest and smallest simplex
			// points and check for termination.
			if (reldif (f[i_max], f[i_min]) < tol)
				{
				// We're done. Put solution in position 0.
				swap_x (0, i_min);
				swap_f (0, i_min);
				return;
				}

			// Check for too many function evaluations.
			if (evalCount >= MAXEVAL)
				{
				throw new TooManyIterationsException
					("MDMinimizationDownhillSimplex.minimize(): Too many function evaluations ("+
					 MAXEVAL+") with no solution");
				}

			// Try reflecting largest point through opposite face of simplex.
			compute_x_trial (-1.0);

			// If reflection gave a new smallest point, try an additional
			// expansion.
			if (f_trial <= f[i_min])
				{
				compute_x_trial (2.0);
				}

			// If reflection was worse than second-largest point, try a
			// one-dimensional contraction of largest point.
			else if (f_trial >= f[i_2ndmax])
				{
				double prev_f_trial = f_trial;
				compute_x_trial (0.5);

				// If one-dimensional contraction gave no improvement, try an
				// all-dimensional contraction toward smallest point.
				if (f_trial >= prev_f_trial)
					{
					double[] x_min = x[i_min];
					for (int i = 0; i <= N; ++ i)
						{
						if (i != i_min)
							{
							double[] x_i = x[i];
							for (int j = 0; j < N; ++ j)
								{
								x_i[j] = 0.5*(x_i[j] + x_min[j]);
								}
							f[i] = fcn.f (x_i);
							}
						}
					evalCount += N;
					compute_x_sum();
					}
				}
			}
		}

// Hidden operations.

	/**
	 * Do debugging. If the <TT>debug</TT> field is true, the
	 * <TT>subclassDebug()</TT> method is called at the beginning of every
	 * iteration. The fields of this object contain the current state of the
	 * algorithm. The fields of this object must not be altered.
	 * <P>
	 * The default implementation of the <TT>subclassDebug()</TT> method does
	 * nothing. A subclass can override the <TT>subclassDebug()</TT> method to
	 * do something, such as print debugging information.
	 *
	 * @param  iter  Iteration number.
	 * @param  eval  Number of function evaluations.
	 */
	protected void subclassDebug
		(int iter,
		 int eval)
		{
		}

	/**
	 * Computes the sum of the simplex points and puts the answer in x_sum.
	 */
	private void compute_x_sum()
		{
		Arrays.fill (x_sum, 0.0);
		for (int i = 0; i <= N; ++ i)
			{
			double[] x_i = x[i];
			for (int j = 0; j < N; ++ j)
				{
				x_sum[j] += x_i[j];
				}
			}
		}

	/**
	 * Returns the relative difference of the arguments.
	 */
	private double reldif
		(double a,
		 double b)
		{
		return 2.0*abs(a - b)/(abs(a) + abs(b) + 1.0e-10);
		}

	/**
	 * Swap the elements of the x array at the given indexes.
	 */
	private void swap_x
		(int i,
		 int j)
		{
		double[] tmp = x[i];
		x[i] = x[j];
		x[j] = tmp;
		}

	/**
	 * Swap the elements of the f array at the given indexes.
	 */
	private void swap_f
		(int i,
		 int j)
		{
		double tmp = f[i];
		f[i] = f[j];
		f[j] = tmp;
		}

	/**
	 * Computes a trial simplex point by extrapolating the largest point across
	 * the opposite face of the simplex by the given factor. Results are stored
	 * in fields x_trial and f_trial. If the trial simplex point's function
	 * value is less than the largest simplex point's function value, replaces
	 * the largest simplex point with the trial simplex point.
	 */
	private void compute_x_trial
		(double factor)
		{
		double a = (1.0 - factor)/N;
		double b = a - factor;
		double[] x_max = x[i_max];
		for (int i = 0; i < N; ++ i)
			{
			x_trial[i] = a*x_sum[i] - b*x_max[i];
			}
		f_trial = fcn.f (x_trial);
		++ evalCount;
		if (f_trial < f[i_max])
			{
			for (int i = 0; i < N; ++ i)
				{
				x_sum[i] = x_sum[i] - x_max[i] + x_trial[i];
				x_max[i] = x_trial[i];
				}
			f[i_max] = f_trial;
			}
		}

// Unit test main program.

//	/**
//	 * Unit test main program.
//	 */
//	public static void main
//		(String[] args)
//		{
//		MDFunction fcn = new MDFunction()
//			{
//			public int argumentLength()
//				{
//				return 2;
//				}
//			public double f (double[] x)
//				{
//				return 2.0 - gaussian (x[0]) - gaussian (x[1]);
//				}
//			private double gaussian (double x)
//				{
//				double d = x - 0.5;
//				return exp (-d*d);
//				}
//			};
//		MDMinimizationDownhillSimplex minimizer =
//			new MDMinimizationDownhillSimplex (fcn)
//				{
//				protected void subclassDebug (int iter, int eval)
//					{
//					System.out.printf ("--------%n");
//					System.out.printf ("Iteration %d, %d evaluations%n",
//						iter, eval);
//					for (int i = 0; i <= N; ++ i)
//						{
//						System.out.printf ("x[%d] =", i);
//						for (int j = 0; j < N; ++ j)
//							{
//							System.out.printf (" %g", x[i][j]);
//							}
//						System.out.printf (", f[%d] = %g%n", i, f[i]);
//						}
//					}
//				};
//		minimizer.debug = true;
//		minimizer.setSimplex (new double[] {0.0, 0.0}, 0.1);
//		minimizer.minimize();
//		System.out.printf ("--------%n");
//		System.out.printf ("Solution: x =");
//		for (int j = 0; j < minimizer.N; ++ j)
//			{
//			System.out.printf (" %g", minimizer.x[0][j]);
//			}
//		System.out.printf (", f = %g%n", minimizer.f[0]);
//		}

	}
