//******************************************************************************
//
// File:    NonLinearLeastSquares.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.NonLinearLeastSquares
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

package edu.rit.numeric;

/**
 * Class NonLinearLeastSquares provides a method for minimizing the sum of the
 * squares of a series of nonlinear functions. There are <I>M</I> functions,
 * each of which has <I>N</I> inputs. These functions are represented by an
 * object that implements interface {@linkplain VectorFunction}. The
 * <TT>solve()</TT> method finds a vector <B>x</B> such that
 * &Sigma;<SUB><I>i</I></SUB>&nbsp;[<I>f</I><SUB><I>i</I></SUB>(<B>x</B>)]<SUP>2</SUP>
 * is minimized. The inputs to and outputs from the <TT>solve()</TT> method are
 * stored in the fields of an instance of class NonLinearLeastSquares. The
 * Levenberg-Marquardt method is used to find the solution.
 * <P>
 * The Java code is a translation of the Fortran subroutine <TT>LMDER</TT> from
 * the MINPACK library. MINPACK was developed by Jorge Mor&eacute;, Burt Garbow,
 * and Ken Hillstrom at Argonne National Laboratory. For further information,
 * see
 * <A HREF="http://www.netlib.org/minpack/" TARGET="_top">http://www.netlib.org/minpack/</A>.
 *
 * @author  Alan Kaminsky
 * @version 10-Jun-2008
 */
public class NonLinearLeastSquares
	{

// Exported data members.

	/**
	 * The nonlinear functions <I>f</I><SUB><I>i</I></SUB>(<B>x</B>) to be
	 * minimized.
	 */
	public final VectorFunction fcn;

	/**
	 * The number of functions.
	 */
	public final int M;

	/**
	 * The number of arguments for each function.
	 */
	public final int N;

	/**
	 * The <I>N</I>-element <B>x</B> vector for the least squares problem. On
	 * input to the <TT>solve()</TT> method, <TT>x</TT> contains an initial
	 * estimate of the solution vector. On output from the <TT>solve()</TT>
	 * method, <TT>x</TT> contains the final solution vector.
	 */
	public final double[] x;

	/**
	 * The <I>M</I>-element result vector. On output from the <TT>solve()</TT>
	 * method, for <I>i</I> = 0 to <I>M</I>&minus;1,
	 * <TT>fvec</TT><SUB><I>i</I></SUB> = <I>f</I><SUB><I>i</I></SUB>(<B>x</B>),
	 * where <B>x</B> is the final solution vector.
	 */
	public final double[] fvec;

	/**
	 * The <I>M</I>&times;<I>N</I>-element Jacobian matrix. On output from the
	 * <TT>solve()</TT> method, the upper <I>N</I>&times;<I>N</I> submatrix of
	 * <TT>fjac</TT> contains an upper triangular matrix <B>R</B> with diagonal
	 * elements of nonincreasing magnitude such that
	 * <P>
	 * <CENTER>
	 * <B>P</B><SUP>T</SUP>&nbsp;(<B>J</B><SUP>T</SUP>&nbsp;<B>J</B>)&nbsp;<B>P</B>&emsp;=&emsp;<B>R</B><SUP>T</SUP>&nbsp;<B>R</B>
	 * </CENTER>
	 * <P>
	 * where <B>P</B> is a permutation matrix and <B>J</B> is the final
	 * calculated Jacobian. Column <I>j</I> of <B>P</B> is column
	 * <TT>ipvt[j]</TT> (see below) of the identity matrix. The lower
	 * trapezoidal part of <TT>fjac</TT> contains information generated during
	 * the computation of <B>R</B>.
	 * <P>
	 * <TT>fjac</TT> is used to calculate the covariance matrix of the solution.
	 * This calculation is not yet implemented.
	 */
	public final double[][] fjac;

	/**
	 * The <I>N</I>-element permutation vector for <TT>fjac</TT>. On output from
	 * the <TT>solve()</TT> method, <TT>ipvt</TT> defines a permutation matrix
	 * <B>P</B> such that <B>JP</B> = <B>QR</B>, where <B>J</B> is the final
	 * calculated Jacobian, <B>Q</B> is orthogonal (not stored), and <B>R</B> is
	 * upper triangular with diagonal elements of nonincreasing magnitude.
	 * Column <I>j</I> of <B>P</B> is column <TT>ipvt[j]</TT> of the identity
	 * matrix.
	 */
	public final int[] ipvt;

	/**
	 * Tolerance. An input to the <TT>solve()</TT> method. Must be &gt;= 0.
	 * Termination occurs when the algorithm estimates either that the relative
	 * error in the sum of squares is at most <TT>tol</TT> or that the relative
	 * error between <TT>x</TT> and the solution is at most <TT>tol</TT>. The
	 * default tolerance is 1&times;10<SUP>&minus;6</SUP>.
	 */
	public double tol = 1.0e-6;

	/**
	 * Information about the outcome. An output of the <TT>solve()</TT> method.
	 * The possible values are:
	 * <UL>
	 * <LI>
	 * 0 -- Improper input parameters. Also, an IllegalArgumentException is
	 * thrown.
	 * <LI>
	 * 1 -- Algorithm estimates that the relative error in the sum of squares is
	 * at most <TT>tol</TT>.
	 * <LI>
	 * 2 -- Algorithm estimates that the relative error between <TT>x</TT> and
	 * the solution is at most <TT>tol</TT>.
	 * <LI>
	 * 3 -- Conditions for <TT>info</TT> = 1 and <TT>info</TT> = 2 both hold.
	 * <LI>
	 * 4 -- <TT>fvec</TT> is orthogonal to the columns of the Jacobian to
	 * machine precision.
	 * <LI>
	 * 5 -- Number of function evaluations has reached 100(<I>N</I>+1). Also, a
	 * TooManyIterationsException is thrown.
	 * <LI>
	 * 6 -- <TT>tol</TT> is too small. No further reduction in the sum of
	 * squares is possible.
	 * <LI>
	 * 7 -- <TT>tol</TT> is too small. No further improvement in the approximate
	 * solution <TT>x</TT> is possible.
	 * </UL>
	 */
	public int info;

	/**
	 * Debug printout flag. An input to the <TT>solve()</TT> method. If
	 * <TT>nprint</TT> &gt; 0, the <TT>subclassDebug()</TT> method is called at
	 * the beginning of the first iteration and every <TT>nprint</TT> iterations
	 * thereafter and immediately prior to return. If <TT>nprint</TT> &lt;= 0,
	 * the <TT>subclassDebug()</TT> method is not called. The default setting is
	 * 0.
	 */
	public int nprint = 0;

	// Working variables.
	private double[] diag;
	private double[] qtf;
	private double[] wa1;
	private double[] wa2;
	private double[] wa3;
	private double[] wa4;
	private int nfev;
	private int njev;

// Hidden constants.

	// Machine epsilon.
	private static final double dpmpar_1 = 2.22044604926e-16;

	// Smallest positive nonzero double value.
	private static final double dpmpar_2 = 2.22507385852e-308;

	// Largest positive double value.
	private static final double dpmpar_3 = 1.79769313485e+308;

	// Used in the enorm() method.
	private static final double rdwarf = 3.834e-20;
	private static final double rgiant = 1.304e+19;

// Exported constructors.

	/**
	 * Construct a new nonlinear least squares problem for the given functions.
	 * Field <TT>fcn</TT> is set to <TT>theFunction</TT>. Fields <TT>M</TT> and
	 * <TT>N</TT> are set by calling the <TT>resultLength()</TT> and
	 * <TT>argumentLength()</TT> methods of <TT>theFunction</TT>. The vector and
	 * matrix fields <TT>x</TT>, <TT>fvec</TT>, <TT>fjac</TT>, and <TT>ipvt</TT>
	 * are allocated with the proper sizes but are not filled in.
	 *
	 * @param  theFunction  Nonlinear functions to be minimized.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFunction</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <I>M</I> &lt;= 0, <I>N</I> &lt;= 0,
	 *     or <I>M</I> &lt; <I>N</I>.
	 */
	public NonLinearLeastSquares
		(VectorFunction theFunction)
		{
		fcn = theFunction;
		M = theFunction.resultLength();
		N = theFunction.argumentLength();
		if (M <= 0)
			{
			throw new IllegalArgumentException
				("NonLinearLeastSquares(): M (= "+M+") <= 0, illegal");
			}
		if (N <= 0)
			{
			throw new IllegalArgumentException
				("NonLinearLeastSquares(): N (= "+N+") <= 0, illegal");
			}
		if (M < N)
			{
			throw new IllegalArgumentException
				("NonLinearLeastSquares(): M (= "+M+") < N (= "+N+
				 "), illegal");
			}

		x = new double [N];
		fvec = new double [M];
		fjac = new double [M] [N];
		ipvt = new int [N];

		diag = new double [N];
		qtf = new double [N];
		wa1 = new double [N];
		wa2 = new double [N];
		wa3 = new double [N];
		wa4 = new double [M];
		}

// Exported operations.

	/**
	 * Solve this nonlinear least squares minimization problem. The
	 * <TT>solve()</TT> method finds a vector <B>x</B> such that
	 * &Sigma;<SUB><I>i</I></SUB>&nbsp;[<I>f</I><SUB><I>i</I></SUB>(<B>x</B>)]<SUP>2</SUP>
	 * is minimized. On input, the field <TT>x</TT> must be filled in with an
	 * initial estimate of the solution vector, and the field <TT>tol</TT> must
	 * be set to the desired tolerance. On output, the other fields are filled
	 * in with the solution as explained in the documentation for each field.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>tol</TT> &lt;= 0.
	 * @exception  TooManyIterationsException
	 *     (unchecked exception) Thrown if too many iterations occurred without
	 *     finding a minimum (100(<I>N</I>+1) iterations).
	 */
	public void solve()
		{
		info = 0;
		if (tol < 0.0)
			{
			throw new IllegalArgumentException
				("NonLinearLeastSquares.solve(): tol (= "+tol+") < 0, illegal");
			}
		lmder
			(/*ftol  */ tol,
			 /*xtol  */ tol,
			 /*gtol  */ 0.0,
			 /*maxfev*/ 100*(N+1),
			 /*mode  */ 1,
			 /*factor*/ 100.0);
		if (info == 5)
			{
			throw new TooManyIterationsException
				("NonLinearLeastSquares.solve(): Too many iterations");
			}
		else if (info == 8)
			{
			info = 4;
			}
		}

// Hidden operations.

	/**
	 * Levenberg-Marquardt algorithm with user-supplied derivatives.
     * <P>
	 * The purpose of lmder() is to minimize the sum of the squares of m
	 * nonlinear functions in n variables by a modification of the
	 * Levenberg-Marquardt algorithm. The user must provide a subroutine which
	 * calculates the functions and the Jacobian.
	 *
	 * @param  ftol
	 *     A nonnegative input variable. Termination occurs when both the actual
	 *     and predicted relative reductions in the sum of squares are at most
	 *     ftol. Therefore, ftol measures the relative error desired in the sum
	 *     of squares.
	 * @param  xtol
	 *     A nonnegative input variable. Termination occurs when the relative
	 *     error between two consecutive iterates is at most xtol. Therefore,
	 *     xtol measures the relative error desired in the approximate solution.
	 * @param  gtol
	 *     A nonnegative input variable. Termination occurs when the cosine of
	 *     the angle between fvec and any column of the jacobian is at most gtol
	 *     in absolute value. Therefore, gtol measures the orthogonality desired
	 *     between the function vector and the columns of the jacobian.
	 * @param  maxfev
	 *     A positive integer input variable. Termination occurs when the number
	 *     of calls to fcn.f() has reached maxfev.
	 * @param  mode
	 *     An integer input variable. If mode = 1, the variables will be scaled
	 *     internally. If mode = 2, the scaling is specified by the input diag.
	 *     Other values of mode are equivalent to mode = 1.
	 * @param  factor
	 *     A positive input variable used in determining the initial step bound.
	 *     This bound is set to the product of factor and the Euclidean norm of
	 *     diag*x if nonzero, or else to factor itself. In most cases factor
	 *     should lie in the interval (0.1,100.0). 100.0 is a generally
	 *     recommended value.
	 */
	private void lmder
		(double ftol,
		 double xtol,
		 double gtol,
		 int maxfev,
		 int mode,
		 double factor)
		{
		int iter, l;
		double actred, delta, dirder, epsmch, fnorm, fnorm1, gnorm, par, pnorm;
		double prered, ratio, sum, temp, temp1, temp2, xnorm;

		// epsmch is the machine precision.
		epsmch = dpmpar_1;

		info = 0;
		nfev = 0;
		njev = 0;
		delta = 0.0;
		xnorm = 0.0;

		// Check the input parameters for errors.
		if (ftol < 0.0)
			{
			throw new IllegalArgumentException
				("NonLinearLeastSquares.lmder(): ftol (= "+ftol+
				 ") < 0, illegal");
			}
		if (xtol < 0.0)
			{
			throw new IllegalArgumentException
				("NonLinearLeastSquares.lmder(): xtol (= "+xtol+
				 ") < 0, illegal");
			}
		if (gtol < 0.0)
			{
			throw new IllegalArgumentException
				("NonLinearLeastSquares.lmder(): gtol (= "+gtol+
				 ") < 0, illegal");
			}
		if (maxfev <= 0)
			{
			throw new IllegalArgumentException
				("NonLinearLeastSquares.lmder(): maxfev (= "+maxfev+
				 ") <= 0, illegal");
			}
		if (factor <= 0.0)
			{
			throw new IllegalArgumentException
				("NonLinearLeastSquares.lmder(): factor (= "+factor+
				 ") <= 0, illegal");
			}
		if (mode == 2)
			{
			for (int j = 0; j < N; ++ j)
				{
				if (diag[j] <= 0.0)
					{
					throw new IllegalArgumentException
						("NonLinearLeastSquares.lmder(): diag["+j+"] (= "+
						 diag[j]+") <= 0, illegal");
					}
				}
			}

		// Evaluate the function at the starting point and calculate its norm.
		fcn.f(x,fvec);
		++ nfev;
		fnorm = enorm(M,fvec);

		// Initialize Levenberg-Marquardt parameter and iteration counter.
		par = 0.0;
		iter = 1;

		// Beginning of the outer loop.
		outerloop: for (;;)
			{
			// Calculate the Jacobian matrix.
			fcn.df(x,fjac);
			++ njev;

			// If requested, print debug information.
			if (nprint > 0 && ((iter-1) % nprint) == 0) subclassDebug (iter);

			// Compute the QR factorization of the Jacobian.
			qrfac(M,N,fjac,true,ipvt,wa1,wa2,wa3);

			if (iter == 1)
				{
				// On the first iteration and if mode is 1, scale according to
				// the norms of the columns of the initial Jacobian.
				if (mode != 2)
					{
					for (int j = 1; j < N; ++ j)
						{
						diag[j] = wa2[j];
						if (wa2[j] == 0.0) diag[j] = 1.0;
						}
					}

				// On the first iteration, calculate the norm of the scaled x
				// and initialize the step bound delta.
				for (int j = 0; j < N; ++ j)
					{
					wa3[j] = diag[j]*x[j];
					}
				xnorm = enorm(N,wa3);
				delta = factor*xnorm;
				if (delta == 0.0) delta = factor;
				}

			// Form (q transpose)*fvec and store the first N components in qtf.
			for (int i = 0; i < M; ++ i)
				{
				wa4[i] = fvec[i];
				}
			for (int j = 0; j < N; ++ j)
				{
				if (fjac[j][j] != 0.0)
					{
					sum = 0.0;
					for (int i = j; i < M; ++ i)
						{
						sum += fjac[i][j]*wa4[i];
						}
					temp = -sum/fjac[j][j];
					for (int i = j; i < M; ++ i)
						{
						wa4[i] += fjac[i][j]*temp;
						}
					}
				fjac[j][j] = wa1[j];
				qtf[j] = wa4[j];
				}

			// Compute the norm of the scaled gradient.
			gnorm = 0.0;
			if (fnorm != 0.0)
				{
				for (int j = 0; j < N; ++ j)
					{
					l = ipvt[j];
					if (wa2[l] != 0.0)
						{
						sum = 0.0;
						for (int i = 0; i <= j; ++ i)
							{
							sum += fjac[i][j]*(qtf[i]/fnorm);
							}
						gnorm = Math.max(gnorm,Math.abs(sum/wa2[l]));
						}
					}
				}

			// Test for convergence of the gradient norm.
			if (gnorm <= gtol)
				{
				info = 4;
				break outerloop;
				}

			// Rescale if necessary.
			if (mode != 2)
				{
				for (int j = 0; j < N; ++ j)
					{
					diag[j] = Math.max(diag[j],wa2[j]);
					}
				}

			// Beginning of the inner loop.
			innerloop: for (;;)
				{
				// Determine the Levenberg-Marquardt parameter.
				par = lmpar(N,fjac,ipvt,diag,qtf,delta,par,wa1,wa2,wa3,wa4);

				// Store the direction p and x + p. Calculate the norm of p.
				for (int j = 0; j < N; ++ j)
					{
					wa1[j] = -wa1[j];
					wa2[j] = x[j] + wa1[j];
					wa3[j] = diag[j]*wa1[j];
					}
				pnorm = enorm(N,wa3);

				// On the first iteration, adjust the initial step bound.
				if (iter == 1) delta = Math.min(delta,pnorm);

				// Evaluate the function at x + p and calculate its norm.
				fcn.f(wa2,wa4);
				++ nfev;
				fnorm1 = enorm(M,wa4);

				// Compute the scaled actual reduction.
				actred = -1.0;
				if (0.1*fnorm1 < fnorm) actred = 1.0 - sqr(fnorm1/fnorm);

				// Compute the scaled predicted reduction and the scaled
				// directional derivative.
				for (int j = 0; j < N; ++ j)
					{
					wa3[j] = 0.0;
					l = ipvt[j];
					temp = wa1[l];
					for (int i = 0; i <= j; ++ i)
						{
						wa3[i] += fjac[i][j]*temp;
						}
					}
				temp1 = enorm(N,wa3)/fnorm;
				temp2 = (Math.sqrt(par)*pnorm)/fnorm;
				prered = sqr(temp1) + sqr(temp2)*2.0;
				dirder = -(sqr(temp1) + sqr(temp2));

				// Compute the ratio of the actual to the predicted reduction.
				ratio = 0.0;
				if (prered != 0.0) ratio = actred/prered;

				// Update the step bound.
				if (ratio <= 0.25)
					{
					temp =
						actred >= 0.0 ?
							0.5 :
							0.5*dirder/(dirder + 0.5*actred);
					if (0.1*fnorm1 >= fnorm || temp < 0.1) temp = 0.1;
					delta = temp*Math.min(delta,pnorm*10.0);
					par = par/temp;
					}
				else if (par == 0.0 || ratio >= 0.75)
					{
					delta = pnorm*2.0;
					par = 0.5*par;
					}

				// Test for successful iteration.
				if (ratio >= 0.0001)
					{
					// Successful iteration. Update x, fvec, and their norms.
					for (int j = 0; j < N; ++ j)
						{
						x[j] = wa2[j];
						wa2[j] = diag[j]*x[j];
						}
					for (int i = 0; i < M; ++ i)
						{
						fvec[i] = wa4[i];
						}
					xnorm = enorm(N,wa2);
					fnorm = fnorm1;
					++ iter;
					}

				// Tests for convergence.
				boolean fconv =
					Math.abs(actred) <= ftol && prered <= ftol && ratio <= 2.0;
				boolean xconv =
					delta <= xtol*xnorm;
				info = (fconv ? 1 : 0) + (xconv ? 2 : 0);
				if (info != 0) break outerloop;

				// Tests for termination and stringent tolerances.
				if (nfev >= maxfev) info = 5;
				if (Math.abs(actred) <= epsmch &&
					prered <= epsmch && ratio <= 2.0) info = 6;
				if (delta <= epsmch*xnorm) info = 7;
				if (gnorm <= epsmch) info = 8;
				if (info != 0) break outerloop;

				// End of the inner loop. Repeat if iteration unsuccessful.
				if (ratio >= 0.0001) break innerloop;
				}

			// End of the outer loop.
			}

		// Finished.
		if (nprint > 0) subclassDebug (iter);
		}

	/**
	 * Calculate the Levenberg-Marquardt parameter.
	 * <P>
	 * Given an M by N matrix a, an N by N nonsingular diagonal matrix d, an
	 * M-vector b, and a positive number delta, the problem is to determine a
	 * value for the parameter par such that if x solves the system
	 * <PRE>
	 *     a*x = b ,     sqrt(par)*d*x = 0 ,
	 * </PRE>
	 * in the least squares sense, and dxnorm is the euclidean norm of d*x, then
	 * either par is zero and
	 * <PRE>
	 *     (dxnorm-delta) .le. 0.1*delta ,
	 * </PRE>
	 * or par is positive and
	 * <PRE>
	 *     abs(dxnorm-delta) .le. 0.1*delta .
	 * </PRE>
	 * <P>
	 * This subroutine completes the solution of the problem if it is provided
	 * with the necessary information from the QR factorization, with column
	 * pivoting, of a. That is, if a*p = q*r, where p is a permutation matrix, q
	 * has orthogonal columns, and r is an upper triangular matrix with diagonal
	 * elements of nonincreasing magnitude, then lmpar() expects the full upper
	 * triangle of r, the permutation matrix p, and the first n components of (q
	 * transpose)*b. On output lmpar() also provides an upper triangular matrix
	 * s such that
	 * <PRE>
	 *      t   t                   t
	 *     p *(a *a + par*d*d)*p = s *s .
	 * </PRE>
	 * s is employed within lmpar and may be of separate interest.
	 * <P>
	 * Only a few iterations are generally needed for convergence of the
	 * algorithm. If, however, the limit of 10 iterations is reached, then the
	 * output par will contain the best value obtained so far.
	 *
	 * @param  n
	 *     A positive integer input variable set to the order of r.
	 * @param  r
	 *     An n by n array. On input the full upper triangle must contain the
	 *     full upper triangle of the matrix r. On output the full upper
	 *     triangle is unaltered, and the strict lower triangle contains the
	 *     strict upper triangle (transposed) of the upper triangular matrix s.
	 * @param  ipvt
	 *     An integer input array of length n which defines the permutation
	 *     matrix p such that a*p = q*r. Column j of p is column ipvt[j] of the
	 *     identity matrix.
	 * @param  diag
	 *     An input array of length n which must contain the diagonal elements
	 *     of the matrix d.
	 * @param  qtb
	 *     An input array of length n which must contain the first n elements of
	 *     the vector (q transpose)*b.
	 * @param  delta
	 *     A positive input variable which specifies an upper bound on the
	 *     Euclidean norm of d*x.
	 * @param  par
	 *     A nonnegative variable. On input par contains an initial estimate of
	 *     the levenberg-marquardt parameter. The return value of lmpar() is the
	 *     final estimate for par.
	 * @param  x
	 *     An output array of length n which contains the least squares solution
	 *     of the system a*x = b, sqrt(par)*d*x = 0, for the returned par.
	 * @param  sdiag
	 *     An output array of length n which contains the diagonal elements of
	 *     the upper triangular matrix s.
	 * @param  wa1
	 *     A work array of length n.
	 * @param  wa2
	 *     A work array of length n.
	 *
	 * @return  Final estimate of the Levenberg-Marquardt parameter par.
	 */
	private static double lmpar
		(int n,
		 double[][] r,
		 int[] ipvt,
		 double[] diag,
		 double[] qtb,
		 double delta,
		 double par,
		 double[] x,
		 double[] sdiag,
		 double[] wa1,
		 double[] wa2)
		{
		int iter, nsing;
		double dxnorm, dwarf, fp, gnorm, parc, parl, paru, sum, temp;

		// dwarf is the smallest positive magnitude.
		dwarf = dpmpar_2;

		// Compute and store in x the Gauss-Newton direction. If the Jacobian is
		// rank-deficient, obtain a least-squares solution.
		nsing = n;
		for (int j = 0; j < n; ++ j)
			{
			wa1[j] = qtb[j];
			if (r[j][j] == 0.0 && nsing == n) nsing = j;
			if (nsing < n) wa1[j] = 0.0;
			}
		for (int j = nsing-1; j >= 0; -- j)
			{
			wa1[j] = wa1[j]/r[j][j];
			temp = wa1[j];
			for (int i = 0; i < j; ++ i)
				{
				wa1[i] -= r[i][j]*temp;
				}
			}
		for (int j = 0; j < n; ++ j)
			{
			x[ipvt[j]] = wa1[j];
			}

		// Initialize the iteration counter. Evaluate the function at the
		// origin, and test for acceptance of the Gauss-Newton direction.
		iter = 0;
		for (int j = 0; j < n; ++ j)
			{
			wa2[j] = diag[j]*x[j];
			}
		dxnorm = enorm(n,wa2);
		fp = dxnorm - delta;
		if (fp <= 0.1*delta) return 0.0;

		// If the Jacobian is not rank deficient, the Newton step provides a
		// lower bound, parl, for the zero of the function. Otherwise set this
		// bound to zero.
		parl = 0.0;
		if (nsing == n)
			{
			for (int j = 0; j < n; ++ j)
				{
				int l = ipvt[j];
				wa1[j] = diag[l]*(wa2[l]/dxnorm);
				}
			for (int j = 0; j < n; ++ j)
				{
				sum = 0.0;
				for (int i = 0; i < j; ++ i)
					{
					sum += r[i][j]*wa1[i];
					}
				wa1[j] = (wa1[j] - sum)/r[j][j];
				}
			temp = enorm(n,wa1);
			parl = ((fp/delta)/temp)/temp;
			}

		// Calculate an upper bound, paru, for the zero of the function.
		for (int j = 0; j < n; ++ j)
			{
			sum = 0.0;
			for (int i = 0; i <= j; ++ i)
				{
				sum += r[i][j]*qtb[i];
				}
			wa1[j] = sum/diag[ipvt[j]];
			}
		gnorm = enorm(n,wa1);
		paru = gnorm/delta;
		if (paru == 0.0) paru = dwarf/Math.min(delta,0.1);

		// If the input par lies outside of the interval (parl,paru), set par to
		// the closer endpoint.
		par = Math.max(par,parl);
		par = Math.min(par,paru);
		if (par == 0.0) par = gnorm/dxnorm;

		iterloop: for (;;)
			{
			// Beginning of an iteration.
			++ iter;

			// Evaluate the function at the current value of par.
			if (par == 0.0) par = Math.max(dwarf,0.001*paru);
			temp = Math.sqrt(par);
			for (int j = 0; j < n; ++ j)
				{
				wa1[j] = temp*diag[j];
				}
			qrsolv(n,r,ipvt,wa1,qtb,x,sdiag,wa2);
			for (int j = 0; j < n; ++ j)
				{
				wa2[j] = diag[j]*x[j];
				}
			dxnorm = enorm(n,wa2);
			temp = fp;
			fp = dxnorm - delta;

			// If the function is small enough, accept the current value of par.
			// Also test for the exceptional cases where parl is zero or the
			// number of iterations has reached 10.
			if (Math.abs(fp) <= 0.1*delta ||
					(parl == 0.0 && fp <= temp && temp < 0.0) ||
					iter == 10)
				{
				break iterloop;
				}

			// Compute the Newton correction.
			for (int j = 0; j < n; ++ j)
				{
				int l = ipvt[j];
				wa1[j] = diag[l]*(wa2[l]/dxnorm);
				}
			for (int j = 0; j < n; ++ j)
				{
				wa1[j] /= sdiag[j];
				temp = wa1[j];
				for (int i = j+1; i < n; ++ i)
					{
					wa1[i] -= r[i][j]*temp;
					}
				}
			temp = enorm(n,wa1);
			parc = ((fp/delta)/temp)/temp;

			// Depending on the sign of the function, update parl or paru.
			if (fp > 0.0) parl = Math.max(parl,par);
			if (fp < 0.0) paru = Math.min(paru,par);

			// Compute an improved estimate for par.
			par = Math.max(parl,par+parc);

			// End of an iteration.
			}

		// Finished.
		return par;
		}

	/**
	 * Compute the QR factorization of a matrix.
	 * <P>
	 * This subroutine uses Householder transformations with column pivoting
	 * (optional) to compute a QR factorization of the m by n matrix a. That is,
	 * qrfac() determines an orthogonal matrix q, a permutation matrix p, and an
	 * upper trapezoidal matrix r with diagonal elements of nonincreasing
	 * magnitude, such that a*p = q*r. The Householder transformation for column
	 * k, k = 1,2,...,min(m,n), is of the form
	 * <PRE>
	 *                     t
	 *     i - (1/u(k))*u*u
	 * </PRE>
	 * where u has zeros in the first k-1 positions. The form of this
	 * transformation and the method of pivoting first appeared in the
	 * corresponding LINPACK subroutine.
	 *
	 * @param  m
	 *     A positive integer input variable set to the number of rows of a.
	 * @param  n
	 *     A positive integer input variable set to the number of columns of a.
	 * @param  a
	 *     An m by n array. On input a contains the matrix for which the QR
	 *     factorization is to be computed. On output the strict upper
	 *     trapezoidal part of a contains the strict upper trapezoidal part of
	 *     r, and the lower trapezoidal part of a contains a factored form of q
	 *     (the non-trivial elements of the u vectors described above).
	 * @param  pivot
	 *     A Boolean input variable. If pivot is set true, then column pivoting
	 *     is enforced. If pivot is set false, then no column pivoting is done.
	 * @param  ipvt
	 *     An integer output array of length n. ipvt defines the permutation
	 *     matrix p such that a*p = q*r. Column j of p is column ipvt(j) of the
	 *     identity matrix. If pivot is false, ipvt is not referenced.
	 * @param  rdiag
	 *     An output array of length n which contains the diagonal elements of
	 *     r.
	 * @param  acnorm
	 *     An output array of length n which contains the norms of the
	 *     corresponding columns of the input matrix a. If this information is
	 *     not needed, then acnorm can coincide with rdiag.
	 * @param  wa
	 *     A work array of length n. If pivot is false, then wa can coincide
	 *     with rdiag.
	 */
	private static void qrfac
		(int m,
		 int n,
		 double[][] a,
		 boolean pivot,
		 int[] ipvt,
		 double[] rdiag,
		 double[] acnorm,
		 double[] wa)
		{
		int kmax, minmn, itemp;
		double ajnorm, epsmch, sum, temp;

		// epsmch is the machine precision.
		epsmch = dpmpar_1;

		// Compute the initial column norms and initialize several arrays.
		for (int j = 0; j < n; ++ j)
			{
			acnorm[j] = enorm(m,a,0,j);
			rdiag[j] = acnorm[j];
			wa[j] = rdiag[j];
			if (pivot) ipvt[j] = j;
			}

		// Reduce a to r with Householder transformations.
		minmn = Math.min(m,n);
		for (int j = 0; j < minmn; ++ j)
			{
			if (pivot)
				{
				// Bring the column of largest norm into the pivot position.
				kmax = j;
				for (int k = j; k < n; ++ k)
					{
					if (rdiag[k] > rdiag[kmax]) kmax = k;
					}
				if (kmax != j)
					{
					for (int i = 0; i < m; ++ i)
						{
						temp = a[i][j];
						a[i][j] = a[i][kmax];
						a[i][kmax] = temp;
						}
					rdiag[kmax] = rdiag[j];
					wa[kmax] = wa[j];
					itemp = ipvt[j];
					ipvt[j] = ipvt[kmax];
					ipvt[kmax] = itemp;
					}
				}

			// Compute the Householder transformation to reduce the j-th column
			// of a to a multiple of the j-th unit vector.
			ajnorm = enorm(m,a,j,j);
			if (ajnorm != 0.0)
				{
				if (a[j][j] < 0.0) ajnorm = -ajnorm;
				for (int i = j; i < m; ++ i)
					{
					a[i][j] /= ajnorm;
					}
				a[j][j] += 1.0;

				// Apply the transformation to the remaining columns and update
				// the norms.
				for (int k = j+1; k < n; ++ k)
					{
					sum = 0.0;
					for (int i = j; i < m; ++ i)
						{
						sum += a[i][j]*a[i][k];
						}
					temp = sum/a[j][j];
					for (int i = j; i < m; ++ i)
						{
						a[i][k] -= temp*a[i][j];
						}
					if (pivot && rdiag[k] != 0.0)
						{
						temp = a[j][k]/rdiag[k];
						rdiag[k] *= Math.sqrt(Math.max(0.0,1.0-sqr(temp)));
						if (0.05*sqr(rdiag[k]/wa[k]) <= epsmch)
							{
							rdiag[k] = enorm(m,a,j+1,k);
							wa[k] = rdiag[k];
							}
						}
					}
				}

			rdiag[j] = -ajnorm;
			}
		}

	/**
	 * Solve a linear system of equations using a QR factorization.
	 * <P>
	 * Given an m by n matrix a, an n by n diagonal matrix d, and an m-vector b,
	 * the problem is to determine an x which solves the system
	 * <PRE>
	 *     a*x = b ,     d*x = 0 ,
	 * </PRE>
	 * in the least squares sense.
	 * <P>
	 * This subroutine completes the solution of the problem if it is provided
	 * with the necessary information from the QR factorization, with column
	 * pivoting, of a. That is, if a*p = q*r, where p is a permutation matrix, q
	 * has orthogonal columns, and r is an upper triangular matrix with diagonal
	 * elements of nonincreasing magnitude, then qrsolv() expects the full upper
	 * triangle of r, the permutation matrix p, and the first n components of (q
	 * transpose)*b. The system a*x = b, d*x = 0, is then equivalent to
	 * <PRE>
	 *            t       t
	 *     r*z = q *b ,  p *d*p*z = 0 ,
	 * </PRE>
	 * where x = p*z. if this system does not have full rank, then a least
	 * squares solution is obtained. On output qrsolv() also provides an upper
	 * triangular matrix s such that
	 * <PRE>
	 *      t   t               t
	 *     p *(a *a + d*d)*p = s *s .
	 * </PRE>
	 * s is computed within qrsolv() and may be of separate interest.
	 *
	 * @param  n
	 *     A positive integer input variable set to the order of r.
	 * @param  r
	 *     An n by n array. On input the full upper triangle must contain the
	 *     full upper triangle of the matrix r. On output the full upper
	 *     triangle is unaltered, and the strict lower triangle contains the
	 *     strict upper triangle (transposed) of the upper triangular matrix s.
	 * @param  ipvt
	 *     An integer input array of length n which defines the permutation
	 *     matrix p such that a*p = q*r. Column j of p is column ipvt(j) of the
	 *     identity matrix.
	 * @param  diag
	 *     An input array of length n which must contain the diagonal elements
	 *     of the matrix d.
	 * @param  qtb
	 *     An input array of length n which must contain the first n elements of
	 *     the vector (q transpose)*b.
	 * @param  x
	 *     An output array of length n which contains the least squares solution
	 *     of the system a*x = b, d*x = 0.
	 * @param  sdiag
	 *     An output array of length n which contains the diagonal elements of
	 *     the upper triangular matrix s.
	 * @param  wa
	 *     A work array of length n.
	 */
	private static void qrsolv
		(int n,
		 double[][] r,
		 int[] ipvt,
		 double[] diag,
		 double[] qtb,
		 double[] x,
		 double[] sdiag,
		 double[] wa)
		{
		int nsing;
		double cos, cotan, qtbpj, sin, sum, tan, temp;

		// Copy r and (q transpose)*b to preserve input and initialize s. In
		// particular, save the diagonal elements of r in x.
		for (int j = 0; j < n; ++ j)
			{
			for (int i = j; i < n; ++ i)
				{
				r[i][j] = r[j][i];
				}
			x[j] = r[j][j];
			wa[j] = qtb[j];
			}

		// Eliminate the diagonal matrix d using a Givens rotation.
		for (int j = 0; j < n; ++ j)
			{
			// Prepare the row of d to be eliminated, locating the diagonal
			// element using p from the QR factorization.
			int l = ipvt[j];
			if (diag[l] != 0.0)
				{
				for (int k = j; k < n; ++ k)
					{
					sdiag[k] = 0.0;
					}
				sdiag[j] = diag[l];

				// The transformations to eliminate the row of d modify only a
				// single element of (q transpose)*b beyond the first n, which
				// is initially zero.
				qtbpj = 0.0;
				for (int k = j; k < n; ++ k)
					{
					// Determine a Givens rotation which eliminates the
					// appropriate element in the current row of d.
					if (sdiag[k] != 0.0)
						{
						if (Math.abs(r[k][k]) < Math.abs(sdiag[k]))
							{
							cotan = r[k][k]/sdiag[k];
							sin = 0.5/Math.sqrt(0.25+0.25*sqr(cotan));
							cos = sin*cotan;
							}
						else
							{
							tan = sdiag[k]/r[k][k];
							cos = 0.5/Math.sqrt(0.25+0.25*sqr(tan));
							sin = cos*tan;
							}

						// Compute the modified diagonal element of r and the
						// modified element of ((q transpose)*b, 0).
						r[k][k] = cos*r[k][k] + sin*sdiag[k];
						temp = cos*wa[k] + sin*qtbpj;
						qtbpj = -sin*wa[k] + cos*qtbpj;
						wa[k] = temp;

						// Accumulate the transformation in the row of s.
						for (int i = k+1; i < n; ++ i)
							{
							temp = cos*r[i][k] + sin*sdiag[i];
							sdiag[i] = -sin*r[i][k] + cos*sdiag[i];
							r[i][k] = temp;
							}
						}
					}
				}

			// Store the diagonal element of s and restore the corresponding
			// diagonal element of r.
			sdiag[j] = r[j][j];
			r[j][j] = x[j];
			}

		// Solve the triangular system for z. If the system is singular, then
		// obtain a least squares solution.
		nsing = n;
		for (int j = 0; j < n; ++ j)
			{
			if (sdiag[j] == 0.0 && nsing == n) nsing = j;
			if (nsing < n) wa[j] = 0.0;
			}
		for (int j = nsing-1; j >= 0; -- j)
			{
			sum = 0.0;
			for (int i = j+1; i < nsing; ++ i)
				{
				sum += r[i][j]*wa[i];
				}
			wa[j] = (wa[j] - sum)/sdiag[j];
			}

		// Permute the components of z back to components of x.
		for (int j = 0; j < n; ++ j)
			{
			x[ipvt[j]] = wa[j];
			}
		}

	/**
	 * Calculate the Euclidean norm of the given vector.
	 *
	 * The Euclidean norm is computed by accumulating the sum of squares in
	 * three different sums. The sums of squares for the small and large
	 * components are scaled so that no overflows occur. Non-destructive
	 * underflows are permitted. Underflows and overflows do not occur in the
	 * computation of the unscaled sum of squares for the intermediate
	 * components. The definitions of small, intermediate and large components
	 * depend on two constants, rdwarf and rgiant. The main restrictions on
	 * these constants are that rdwarf**2 not underflow and rgiant**2 not
	 * overflow. The constants given here are suitable for every known computer.
	 *
	 * @param  n  Number of elements in the vector.
	 * @param  x  Vector.
	 *
	 * @return  Euclidean norm of <TT>x</TT>.
	 */
	private static double enorm
		(int n,
		 double[] x)
		{
		double s1 = 0.0;
		double s2 = 0.0;
		double s3 = 0.0;
		double x1max = 0.0;
		double x3max = 0.0;
		double agiant = rgiant/n;
		for (int i = 0; i < n; ++ i)
			{
			double xabs = Math.abs (x[i]);

			// Sum for large components.
			if (xabs >= agiant)
				{
				if (xabs > x1max)
					{
					s1 = 1.0 + s1*sqr(x1max/xabs);
					x1max = xabs;
					}
				else
					{
					s1 += sqr(xabs/x1max);
					}
				}

			// Sum for small components.
			else if (xabs <= rdwarf)
				{
				if (xabs > x3max)
					{
					s3 = 1.0 + s3*sqr(x3max/xabs);
					x3max = xabs;
					}
				else
					{
					if (xabs != 0.0) s3 += sqr(xabs/x3max);
					}
				}

			// Sum for intermediate components.
			else
				{
				s2 += sqr(xabs);
				}
			}

		// Calculation of norm.
		if (s1 != 0.0)
			{
			return x1max*Math.sqrt(s1+(s2/x1max)/x1max);
			}
		else if (s2 != 0.0)
			{
			if (s2 >= x3max)
				{
				return Math.sqrt(s2*(1.0+(x3max/s2)*(x3max*s3)));
				}
			else
				{
				return Math.sqrt(x3max*((s2/x3max)+(x3max*s3)));
				}
			}
		else
			{
			return x3max*Math.sqrt(s3);
			}
		}

	/**
	 * Calculate the Euclidean norm of a portion of one column of the given
	 * matrix.
	 *
	 * @param  n  Number of rows in the matrix.
	 * @param  x  Matrix.
	 * @param  r  Initial row index.
	 * @param  j  Column index.
	 *
	 * @return  Euclidean norm of rows <TT>r</TT> through <TT>n</TT>-1 of column
	 *          <TT>j</TT> of matrix <TT>x</TT>.
	 */
	private static double enorm
		(int n,
		 double[][] x,
		 int r,
		 int j)
		{
		double s1 = 0.0;
		double s2 = 0.0;
		double s3 = 0.0;
		double x1max = 0.0;
		double x3max = 0.0;
		double agiant = rgiant/n;
		for (int i = r; i < n; ++ i)
			{
			double xabs = Math.abs (x[i][j]);

			// Sum for large components.
			if (xabs >= agiant)
				{
				if (xabs > x1max)
					{
					s1 = 1.0 + s1*sqr(x1max/xabs);
					x1max = xabs;
					}
				else
					{
					s1 += sqr(xabs/x1max);
					}
				}

			// Sum for small components.
			else if (xabs <= rdwarf)
				{
				if (xabs > x3max)
					{
					s3 = 1.0 + s3*sqr(x3max/xabs);
					x3max = xabs;
					}
				else
					{
					if (xabs != 0.0) s3 += sqr(xabs/x3max);
					}
				}

			// Sum for intermediate components.
			else
				{
				s2 += sqr(xabs);
				}
			}

		// Calculation of norm.
		if (s1 != 0.0)
			{
			return x1max*Math.sqrt(s1+(s2/x1max)/x1max);
			}
		else if (s2 != 0.0)
			{
			if (s2 >= x3max)
				{
				return Math.sqrt(s2*(1.0+(x3max/s2)*(x3max*s3)));
				}
			else
				{
				return Math.sqrt(x3max*((s2/x3max)+(x3max*s3)));
				}
			}
		else
			{
			return x3max*Math.sqrt(s3);
			}
		}

	/**
	 * Returns the square of x.
	 */
	private static double sqr
		(double x)
		{
		return x*x;
		}

	/**
	 * Print debugging information. If the <TT>nprint</TT> field is greater than
	 * zero, the <TT>subclassDebug()</TT> method is called at the beginning of
	 * the first iteration and every <TT>nprint</TT> iterations thereafter and
	 * immediately prior to return. The fields of this object contain the
	 * current state of the algorithm. The fields of this object must not be
	 * altered.
	 * <P>
	 * The default implementation of the <TT>subclassDebug()</TT> method does
	 * nothing. A subclass can override the <TT>subclassDebug()</TT> method to
	 * do something, such as print debugging information.
	 *
	 * @param  iter  Iteration number.
	 */
	protected void subclassDebug
		(int iter)
		{
		}

	}
