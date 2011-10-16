//******************************************************************************
//
// File:    LinearSolve.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.LinearSolve
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

import java.util.Arrays;

/**
 * Class LinearSolve provides an object for solving a system of linear equations
 * using LU decomposition. The <TT>solve()</TT> method finds a solution to the
 * system of linear equations <B>Ax</B> = <B>b</B>, where <B>A</B> is a square
 * matrix supplied to the constructor and <B>b</B> is a vector supplied to the
 * <TT>solve()</TT> method. Thus, an instance of class LinearSolve can be used
 * to solve many linear systems with the same left-hand-side matrix and
 * different right-hand-side vectors.
 * <P>
 * The Java code for LU decomposition was translated from routine
 * <TT>gsl_linalg_LU_decomp()</TT> in the GNU Scientific Library.
 *
 * @author  Alan Kaminsky
 * @version 07-Jul-2007
 */
public class LinearSolve
	{

// Exported data members.

	/**
	 * The number of rows and columns in the matrix.
	 */
	private int N;

	/**
	 * The <I>N</I>x<I>N</I>-element left-hand-side matrix, containing the LU
	 * decomposition of a rowwise permutation of the original <B>A</B> matrix.
	 */
	private double[][] LU;

	/**
	 * The <I>N</I>-element permutation vector.
	 */
	private int[] p;

	/**
	 * The sign of the permutation.
	 */
	private double signum;

// Exported constructors.

	/**
	 * Construct a new LinearSolve object. <TT>A</TT> must be an
	 * <I>N</I>-by-<I>N</I> matrix with <I>N</I> &gt; 0. This constructor
	 * calculates the LU decomposition of <TT>A</TT> and stores the result
	 * internally; <TT>A</TT> is unchanged.
	 *
	 * @param  A  Left-hand-side matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>A</TT> or any row thereof is
	 *     null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>A</TT> is zero length or is not a
	 *     square matrix.
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if a zero pivot was encountered during
	 *     the LU decomposition.
	 */
	public LinearSolve
		(double[][] A)
		{
		// Copy matrix A and verify preconditions.
		N = A.length;
		if (N == 0)
			{
			throw new IllegalArgumentException
				("LinearSolve(): A is zero length");
			}
		LU = new double [N] [N];
		for (int i = 0; i < N; ++ i)
			{
			if (A[i].length != N)
				{
				throw new IllegalArgumentException
					("LinearSolve(): A is not a square matrix");
				}
			System.arraycopy (A[i], 0, LU[i], 0, N);
			}

		// Allocate storage for permutation.
		p = new int [N];

		// Do the LU decomposition.
		signum = luDecompose (LU, p, N);
		}

// Exported operations.

	/**
	 * Solve the linear system <B>Ax</B> = <B>b</B>. <B>A</B> is the
	 * <I>N</I>-by-<I>N</I> matrix supplied to the constructor. <B>b</B> must be
	 * an <I>N</I>-element array initialized to the right-hand-side vector. The
	 * solution vector is stored in the <I>N</I>-element array <B>x</B>.
	 * <B>b</B> and <B>x</B> must be different arrays.
	 *
	 * @param  x  Solution vector (output).
	 * @param  b  Right-hand-side vector (input).
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>x</TT> or <TT>b</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>x</TT>'s length or <TT>b</TT>'s
	 *     length is not <I>N</I>.
	 */
	public void solve
		(double[] x,
		 double[] b)
		{
		// Verify preconditions.
		if (x.length != N)
			{
			throw new IllegalArgumentException
				("LinearSolve.solve(): x is not a " + N + "-element array");
			}
		if (b.length != N)
			{
			throw new IllegalArgumentException
				("LinearSolve.solve(): b is not a " + N + "-element array");
			}

		// Apply the permutation to b.
		for (int i = 0; i < N; ++ i)
			{
			x[i] = b[p[i]];
			}

		// Compute the solution.
		luSolve (LU, p, x, N);
		}

	/**
	 * Compute <B>A</B><SUP>-1</SUP>, the inverse of <B>A</B>. <B>A</B> is the
	 * <I>N</I>-by-<I>N</I> matrix supplied to the constructor. The inverse of
	 * <B>A</B> is stored in the <I>N</I>-by-<I>N</I> matrix <TT>Ainv</TT>.
	 *
	 * @param  Ainv  Inverse matrix (output).
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>Ainv</TT> or any row thereof is
	 *     null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>Ainv</TT> is not an
	 *     <I>N</I>-by-<I>N</I> matrix.
	 */
	public void invert
		(double[][] Ainv)
		{
		// Verify preconditions.
		if (Ainv.length != N)
			{
			throw new IllegalArgumentException
				("LinearSolve.invert(): Ainv is not a " + N + "-by-" + N +
				 " matrix");
			}
		for (int i = 0; i < N; ++ i)
			{
			if (Ainv[i].length != N)
				{
				throw new IllegalArgumentException
					("LinearSolve.invert(): Ainv is not a " + N + "-by-" + N +
					 " matrix");
				}
			}

		// Allocate temporary storage.
		double[] x = new double [N];

		// Compute and store the columns of the inverse.
		for (int i = 0; i < N; ++ i)
			{
			// Apply the permutation to column i of an identity matrix.
			for (int j = 0; j < N; ++ j)
				{
				x[j] = i == p[j] ? 1.0 : 0.0;
				}
			luSolve (LU, p, x, N);
			for (int j = 0; j < N; ++ j)
				{
				Ainv[j][i] = x[j];
				}
			}
		}

	/**
	 * Compute <B>A</B><SUP>-1</SUP><B>B</B>. <B>A</B> is the
	 * <I>N</I>-by-<I>N</I> matrix supplied to the constructor. <B>B</B> must be
	 * an <I>N</I>-by-<I>N</I> matrix. The matrix product of
	 * <B>A</B><SUP>-1</SUP> and <B>B</B> is stored in the <I>N</I>-by-<I>N</I>
	 * matrix <TT>AinvB</TT>.
	 *
	 * @param  AinvB  Product matrix (output).
	 * @param  B      Matrix (input).
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>AinvB</TT>, <TT>B</TT>, or any
	 *     row thereof is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>AinvB</TT> or <TT>B</TT> is not
	 *     an <I>N</I>-by-<I>N</I> matrix.
	 */
	public void invertMultiply
		(double[][] AinvB,
		 double[][] B)
		{
		// Verify preconditions.
		if (AinvB.length != N)
			{
			throw new IllegalArgumentException
				("LinearSolve.invertMultiply(): AinvB is not a " +
				 N + "-by-" + N + " matrix");
			}
		if (B.length != N)
			{
			throw new IllegalArgumentException
				("LinearSolve.invertMultiply(): B is not a " +
				 N + "-by-" + N + " matrix");
			}
		for (int i = 0; i < N; ++ i)
			{
			if (AinvB[i].length != N)
				{
				throw new IllegalArgumentException
					("LinearSolve.invertMultiply(): AinvB is not a " +
					 N + "-by-" + N + " matrix");
				}
			if (B[i].length != N)
				{
				throw new IllegalArgumentException
					("LinearSolve.invertMultiply(): B is not a " +
					 N + "-by-" + N + " matrix");
				}
			}

		// Allocate temporary storage.
		double[] x = new double [N];

		// Compute and store the columns of the inverse product.
		for (int i = 0; i < N; ++ i)
			{
			// Apply the permutation to column i of B.
			for (int j = 0; j < N; ++ j)
				{
				x[j] = B[p[j]][i];
				}
			luSolve (LU, p, x, N);
			for (int j = 0; j < N; ++ j)
				{
				AinvB[j][i] = x[j];
				}
			}
		}

	/**
	 * Compute det <B>A</B>, the determinant of <B>A</B>. <B>A</B> is the
	 * <I>N</I>-by-<I>N</I> matrix supplied to the constructor. Note that for
	 * larger matrices, det <B>A</B> may overflow or underflow the dynamic range
	 * of type <TT>double</TT>.
	 *
	 * @return  The determinant of <B>A</B>.
	 */
	public double determinant()
		{
		double det = signum;
		for (int i = 0; i < N; ++ i)
			{
			det *= LU[i][i];
			}
		return det;
		}

// Hidden operations.

	/**
	 * Calculate the LU decomposition of matrix A. On input, A must be an NxN
	 * matrix, and P must be an N-element array. On output, A has been replaced
	 * with the LU decomposition of A, and P has been replaced by a description
	 * of a row permutation of A. The upper triangular factor, U, replaces the
	 * diagonal and upper triangle of A. The lower triangular factor, L,
	 * replaces the lower triangle of A; the diagonal elements of L are all 1
	 * and are not stored. The return value is -1 if the row permutation has an
	 * odd number of interchanges or +1 if the row permutation has an even
	 * number of interchanges. The Java code was translated from routine
	 * <TT>gsl_linalg_LU_decomp()</TT> in the GNU Scientific Library.
	 *
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if a zero pivot was encountered during
	 *     the LU decomposition.
	 */
	private static double luDecompose
		(double[][] A,
		 int[] p,
		 int N)
		{
		// Initialize sign of permutation.
		double signum = 1.0;

		// Initialize permutation.
		for (int i = 0; i < N; ++ i)
			{
			p[i] = i;
			}

		// Do all columns.
		for (int j = 0; j < N-1; ++ j)
			{
			// Find pivot element (maximum element) in the j-th column.
			double max = Math.abs (A[j][j]);
			int i_pivot = j;
			for (int i = j+1; i < N; ++ i)
				{
				double aij = Math.abs (A[i][j]);
				if (aij > max)
					{
					max = aij;
					i_pivot = i;
					}
				}

			// If the pivot element is not on the diagonal, interchange rows.
			if (i_pivot != j)
				{
				// Swap pivot row with diagonal row.
				double[] swap = A[i_pivot];
				A[i_pivot] = A[j];
				A[j] = swap;

				// Update permutation.
				int swap2 = p[i_pivot];
				p[i_pivot] = p[j];
				p[j] = swap2;
				signum = - signum;
				}

			// Update the decomposition.
			double ajj = A[j][j];
			if (ajj != 0.0)
				{
				for (int i = j+1; i < N; ++ i)
					{
					double aij = A[i][j] / ajj;
					A[i][j] = aij;
					for (int k = j+1; k < N; ++ k)
						{
						A[i][k] -= aij * A[j][k];
						}
					}
				}
			else // (ajj == 0.0)
				{
				throw new DomainException
					("LinearSolve(): Zero pivot encountered");
				}
			}

		// Return sign of permutation.
		return signum;
		}

	/**
	 * Solve the linear system Ax = b using LU decomposition. On input, LU must
	 * be an NxN matrix which is the output of luDecompose(A,p,N), p must be an
	 * N-element array which is the output of luDecompose(A,p,N), and x must be
	 * an N-element array which is initialized to the right-hand side vector
	 * permuted according to the permutation vector p. On output, LU, p, and b
	 * are unchanged, and x has been replaced with the solution vector.
	 */
	private static void luSolve
		(double[][] LU,
		 int[] p,
		 double[] x,
		 int N)
		{
		// Solve Ly = b using forward substitution. (y uses the same storage as
		// x.)
		for (int i = 1; i < N; ++ i)
			{
			double sum = x[i];
			for (int j = 0; j < i; ++ j)
				{
				sum -= LU[i][j] * x[j];
				}
			x[i] = sum;
			}

		// Solve Ux = y using back substitution.
		x[N-1] /= LU[N-1][N-1];
		for (int i = N-2; i >= 0; -- i)
			{
			double sum = x[i];
			for (int j = i+1; j < N; ++ j)
				{
				sum -= LU[i][j] * x[j];
				}
			x[i] = sum / LU[i][i];
			}
		}

// Unit test main program.

//	/**
//	 * Unit test main program.
//	 */
//	public static void main
//		(String[] args)
//		{
//		int N = Integer.parseInt (args[0]);
//		double[][] A = new double [N] [N];
//		double x = 1.0;
//		System.out.println ("A =");
//		for (int i = 0; i < N; ++ i)
//			{
//			for (int j = 0; j < N; ++ j)
//				{
//				A[i][j] = (i + j) % N + 1;
//				System.out.print (A[i][j]);
//				System.out.print ('\t');
//				x += 1.0;
//				}
//			System.out.println();
//			}
//		LinearSolve solver = new LinearSolve (A);
//		solver.invert (A);
//		System.out.println ("A^{-1} =");
//		for (int i = 0; i < N; ++ i)
//			{
//			for (int j = 0; j < N; ++ j)
//				{
//				System.out.print (A[i][j]);
//				System.out.print ('\t');
//				}
//			System.out.println();
//			}
//		}

	}
