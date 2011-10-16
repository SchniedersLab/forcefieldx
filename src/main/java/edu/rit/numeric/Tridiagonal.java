//******************************************************************************
//
// File:    Tridiagonal.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.Tridiagonal
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

package edu.rit.numeric;

/**
 * Class Tridiagonal provides static methods for solving tridiagonal systems of
 * linear equations.
 * <P>
 * The Java code was translated from routines
 * <TT>gsl_linalg_solve_tridiag()</TT>,
 * <TT>gsl_linalg_solve_symm_tridiag()</TT>,
 * <TT>gsl_linalg_solve_cyc_tridiag()</TT>, and
 * <TT>gsl_linalg_solve_symm_cyc_tridiag()</TT>
 * in the GNU Scientific Library Version 1.9.
 *
 * @author  Alan Kaminsky
 * @version 07-Jul-2007
 */
public class Tridiagonal
	{

// Prevent construction.

	private Tridiagonal()
		{
		}

// Exported operations.

	/**
	 * Solve the given tridiagonal system of linear equations. This method
	 * solves the general <I>N</I>-by-<I>N</I> system <I>Ax</I> = <I>b</I> where
	 * <I>A</I> is tridiagonal (<I>N</I> &gt;= 2). The form of <I>A</I> for the
	 * 4-by-4 case is:
	 * <PRE>
	 *     [ d0  e0  0   0  ]
	 * A = [ f0  d1  e1  0  ]
	 *     [ 0   f1  d2  e2 ]
	 *     [ 0   0   f2  d3 ]
	 * </PRE>
	 *
	 * @param  d  (input) Vector of diagonal elements. Length <I>N</I> must be
	 *            &gt;= 2.
	 * @param  e  (input) Vector of super-diagonal elements. Length must be
	 *            <I>N</I>-1.
	 * @param  f  (input) Vector of sub-diagonal elements. Length must be
	 *            <I>N</I>-1.
	 * @param  b  (input) Vector of right hand side elements. Length must be
	 *            <I>N</I>.
	 * @param  x  (output) Solution vector. Length must be <I>N</I>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if any argument is the wrong length.
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the linear system cannot be solved.
	 */
	public static void solve
		(double[] d,
		 double[] e,
		 double[] f,
		 double[] b,
		 double[] x)
		{
		// Verify preconditions.
		int N = d.length;
		if (N < 2)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solve(): d.length = " + d.length + " illegal");
			}
		if (e.length != N-1)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solve(): e.length = " + e.length + " illegal");
			}
		if (f.length != N-1)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solve(): f.length = " + f.length + " illegal");
			}
		if (b.length != N)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solve(): b.length = " + b.length + " illegal");
			}
		if (x.length != N)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solve(): x.length = " + x.length + " illegal");
			}

		// Working storage.
		double[] alpha = new double [N];
		double[] z = new double [N];

		// Elimination of sub-diagonal. alpha = new diagonal, z = new right hand
		// side.
		alpha[0] = d[0];
		z[0] = b[0];
		if (alpha[0] == 0.0)
			{
			throw new DomainException
				("Tridiagonal.solve(): Zero on diagonal");
			}
		for (int i = 1; i < N; ++ i)
			{
			double t = f[i-1] / alpha[i-1];
			alpha[i] = d[i] - t * e[i-1];
			z[i] = b[i] - t * z[i-1];
			if (alpha[i] == 0.0)
				{
				throw new DomainException
					("Tridiagonal.solve(): Zero on diagonal");
				}
			}

		// Back substitution.
		int Nminus1 = N - 1;
		x[Nminus1] = z[Nminus1] / alpha[Nminus1];
		for (int i = N-2; i >= 0; -- i)
			{
			x[i] = (z[i] - e[i] * x[i+1]) / alpha[i];
			}
		}

	/**
	 * Solve the given symmetric tridiagonal system of linear equations. This
	 * method solves the general <I>N</I>-by-<I>N</I> system <I>Ax</I> =
	 * <I>b</I> where <I>A</I> is symmetric tridiagonal (<I>N</I> &gt;= 2). The
	 * form of <I>A</I> for the 4-by-4 case is:
	 * <PRE>
	 *     [ d0  e0  0   0  ]
	 * A = [ e0  d1  e1  0  ]
	 *     [ 0   e1  d2  e2 ]
	 *     [ 0   0   e2  d3 ]
	 * </PRE>
	 *
	 * @param  d  (input) Vector of diagonal elements. Length <I>N</I> must be
	 *            &gt;= 2.
	 * @param  e  (input) Vector of off-diagonal elements. Length must be
	 *            <I>N</I>-1.
	 * @param  b  (input) Vector of right hand side elements. Length must be
	 *            <I>N</I>.
	 * @param  x  (output) Solution vector. Length must be <I>N</I>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if any argument is the wrong length.
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the linear system cannot be solved.
	 */
	public static void solveSymmetric
		(double[] d,
		 double[] e,
		 double[] b,
		 double[] x)
		{
		solve (d, e, e, b, x);
		}

	/**
	 * Solve the given cyclic tridiagonal system of linear equations. This
	 * method solves the general <I>N</I>-by-<I>N</I> system <I>Ax</I> =
	 * <I>b</I> where <I>A</I> is cyclic tridiagonal (<I>N</I> &gt;= 3). The
	 * form of <I>A</I> for the 4-by-4 case is:
	 * <PRE>
	 *     [ d0  e0  0   f3 ]
	 * A = [ f0  d1  e1  0  ]
	 *     [ 0   f1  d2  e2 ]
	 *     [ e3  0   f2  d3 ]
	 * </PRE>
	 *
	 * @param  d  (input) Vector of diagonal elements. Length <I>N</I> must be
	 *            &gt;= 3.
	 * @param  e  (input) Vector of super-diagonal elements. Length must be
	 *            <I>N</I>.
	 * @param  f  (input) Vector of sub-diagonal elements. Length must be
	 *            <I>N</I>.
	 * @param  b  (input) Vector of right hand side elements. Length must be
	 *            <I>N</I>.
	 * @param  x  (output) Solution vector. Length must be <I>N</I>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if any argument is the wrong length.
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the linear system cannot be solved.
	 */
	public static void solveCyclic
		(double[] d,
		 double[] e,
		 double[] f,
		 double[] b,
		 double[] x)
		{
		// Verify preconditions.
		int N = d.length;
		if (N < 3)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solveCyclic(): d.length = " + d.length +
				 " illegal");
			}
		if (e.length != N)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solveCyclic(): e.length = " + e.length +
				 " illegal");
			}
		if (f.length != N)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solveCyclic(): f.length = " + f.length +
				 " illegal");
			}
		if (b.length != N)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solveCyclic(): b.length = " + b.length +
				 " illegal");
			}
		if (x.length != N)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solveCyclic(): x.length = " + x.length +
				 " illegal");
			}

		// Working storage.
		double[] alpha = new double [N];
		double[] zb = new double [N];
		double[] zu = new double [N];
		double[] w = new double [N];
		double beta;

		// Elimination of sub-diagonal. alpha = new diagonal, zb = new right
		// hand side. A*q = zu.
		if (d[0] == 0.0 || d[1] == 0.0)
			{
			throw new DomainException
				("Tridiagonal.solveCyclic(): Zero on diagonal");
			}
		zb[0] = b[0];
		beta = -d[0];
		double q = 1.0 - (e[0] * f[0]) / (d[0] * d[1]);
		double abs_q_over_beta = Math.abs (q / beta);
		if (abs_q_over_beta <= 0.5)
			{
			}
		else if (abs_q_over_beta < 1.0)
			{
			beta *= 0.5;
			}
		else if (abs_q_over_beta < 2.0)
			{
			beta *= 2.0;
			}
		zu[0] = beta;
		alpha[0] = d[0] - beta;
		if (alpha[0] == 0.0)
			{
			throw new DomainException
				("Tridiagonal.solveCyclic(): Zero on diagonal");
			}
		int Nminus1 = N - 1;
		for (int i = 1; i < Nminus1; ++ i)
			{
			double t = f[i-1] / alpha[i-1];
			alpha[i] = d[i] - t * e[i-1];
			zb[i] = b[i] - t * zb[i-1];
			zu[i] = -t * zu[i-1];
			if (alpha[i] == 0.0)
				{
				throw new DomainException
					("Tridiagonal.solveCyclic(): Zero on diagonal");
				}
			}
		int Nminus2 = N - 2;
		double t = f[Nminus2] / alpha[Nminus2];
		alpha[Nminus1] =
			d[Nminus1] - e[Nminus1] * f[Nminus1] / beta - t * e[Nminus2];
		zb[Nminus1] = b[Nminus1] - t * zb[Nminus2];
		zu[Nminus1] = e[Nminus1] - t * zu[Nminus2];
		if (alpha[Nminus1] == 0.0)
			{
			throw new DomainException
				("Tridiagonal.solveCyclic(): Zero on diagonal");
			}

		// Back substitution.
		w[Nminus1] = zu[Nminus1] / alpha[Nminus1];
		x[Nminus1] = zb[Nminus1] / alpha[Nminus1];
		for (int i = Nminus2; i >= 0; -- i)
			{
			w[i] = (zu[i] - e[i] * w[i+1]) / alpha[i];
			x[i] = (zb[i] - e[i] * x[i+1]) / alpha[i];
			}

		// Sherman-Morrison to fix up from corner elements.
		double vw = w[0] + f[Nminus1] / beta * w[Nminus1] + 1.0;
		double vx = x[0] + f[Nminus1] / beta * x[Nminus1];
		if (vw == 0.0)
			{
			throw new DomainException
				("Tridiagonal.solveCyclic(): Zero on diagonal");
			}
		double vx_over_vw = vx / vw;
		for (int i = 0; i < N; ++ i)
			{
			x[i] -= vx_over_vw * w[i];
			}
		}

	/**
	 * Solve the given symmetric cyclic tridiagonal system of linear equations.
	 * This method solves the general <I>N</I>-by-<I>N</I> system <I>Ax</I> =
	 * <I>b</I> where <I>A</I> is symmetric cyclic tridiagonal (<I>N</I> &gt;=
	 * 3). The form of <I>A</I> for the 4-by-4 case is:
	 * <PRE>
	 *     [ d0  e0  0   e3 ]
	 * A = [ e0  d1  e1  0  ]
	 *     [ 0   e1  d2  e2 ]
	 *     [ e3  0   e2  d3 ]
	 * </PRE>
	 *
	 * @param  d  (input) Vector of diagonal elements. Length <I>N</I> must be
	 *            &gt;= 3.
	 * @param  e  (input) Vector of off-diagonal elements. Length must be
	 *            <I>N</I>.
	 * @param  b  (input) Vector of right hand side elements. Length must be
	 *            <I>N</I>.
	 * @param  x  (output) Solution vector. Length must be <I>N</I>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if any argument is the wrong length.
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the linear system cannot be solved.
	 */
	public static void solveSymmetricCyclic
		(double[] d,
		 double[] e,
		 double[] b,
		 double[] x)
		{
		// Verify preconditions.
		int N = d.length;
		if (N < 3)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solveSymmetricCyclic(): d.length = " + d.length +
				 " illegal");
			}
		if (e.length != N)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solveSymmetricCyclic(): e.length = " + e.length +
				 " illegal");
			}
		if (b.length != N)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solveSymmetricCyclic(): b.length = " + b.length +
				 " illegal");
			}
		if (x.length != N)
			{
			throw new IllegalArgumentException
				("Tridiagonal.solveSymmetricCyclic(): x.length = " + x.length +
				 " illegal");
			}

		// Working storage.
		double[] alpha = new double [N];
		double[] gamma = new double [N];
		double[] delta = new double [N];
		double[] c = new double [N];
		double[] z = new double [N];
		double sum = 0.0;

		// Factor.
		int Nminus1 = N - 1;
		int Nminus2 = N - 2;
		int Nminus3 = N - 3;
		if (d[0] == 0.0)
			{
			throw new DomainException
				("Tridiagonal.solveSymmetricCyclic(): Zero on diagonal");
			}
		alpha[0] = d[0];
		gamma[0] = e[0] / alpha[0];
		delta[0] = e[Nminus1] / alpha[0];
		sum += alpha[0] * delta[0] * delta[0];
		for (int i = 1; i < Nminus2; ++ i)
			{
			alpha[i] = d[i] - e[i-1] * gamma[i-1];
			if (alpha[i] == 0.0)
				{
				throw new DomainException
					("Tridiagonal.solveSymmetricCyclic(): Zero on diagonal");
				}
			gamma[i] = e[i] / alpha[i];
			delta[i] = -delta[i-1] * e[i-1] / alpha[i];
			sum += alpha[i] * delta[i] * delta[i];
			}
		alpha[Nminus2] =
			d[Nminus2] - e[Nminus3] * gamma[Nminus3];
		gamma[Nminus2] =
			(e[Nminus2] - e[Nminus3] * delta[Nminus3]) / alpha[Nminus2];
		alpha[Nminus1] =
			d[Nminus1] - sum - alpha[Nminus2] * gamma[Nminus2] * gamma[Nminus2];

		// Update.
		z[0] = b[0];
		for (int i = 1; i < Nminus1; ++ i)
			{
			z[i] = b[i] - z[i-1] * gamma[i-1];
			}
		sum = 0.0;
		for (int i = 0; i < Nminus2; ++ i)
			{
			sum += delta[i] * z[i];
			}
		z[Nminus1] = b[Nminus1] - sum - gamma[Nminus2] * z[Nminus2];
		for (int i = 0; i < N; ++ i)
			{
			c[i] = z[i] / alpha[i];
			}

		// Back substitution.
		x[Nminus1] = c[Nminus1];
		x[Nminus2] = c[Nminus2] - gamma[Nminus2] * x[Nminus1];
		for (int i = Nminus3; i >= 0; -- i)
			{
			x[i] = c[i] - gamma[i] * x[i+1] - delta[i] * x[Nminus1];
			}
		}

	}
