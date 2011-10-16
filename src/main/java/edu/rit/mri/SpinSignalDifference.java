//******************************************************************************
//
// File:    SpinSignalDifference.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.SpinSignalDifference
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

import edu.rit.numeric.Series;
import edu.rit.numeric.VectorFunction;
import edu.rit.numeric.XYSeries;

/**
 * Class SpinSignalDifference provides an object that computes the difference
 * between a measured spin signal and a model spin signal. Class
 * SpinSignalDifference implements interface {@linkplain
 * edu.rit.numeric.VectorFunction VectorFunction}.
 * <P>
 * An instance of class SpinSignalDifference is constructed with two data
 * series, each an instance of class {@linkplain edu.rit.numeric.Series Series}.
 * The first data series contains the measurement times
 * <I>t</I><SUB><I>i</I></SUB>. The second data series contains the measured
 * spin signals <I>S</I>(<I>t</I><SUB><I>i</I></SUB>). The vector function's
 * result vector length, <I>M</I>, is the same as the length of the data series.
 * <P>
 * The vector function's argument vector <B>x</B> gives parameters for the model
 * spin signal. The argument vector consists of one or more pairs of consecutive
 * values. Each pair of consecutive values corresponds to one tissue. The first
 * value of the pair is <I>&rho;</I><SUB><I>j</I></SUB>, the spin density for
 * tissue <I>j</I>. The second value of the pair is <I>R</I><SUB><I>j</I></SUB>,
 * the spin relaxation rate for tissue <I>j</I>. The number of tissues is
 * specified as a constructor parameter. The vector function's argument vector
 * length, <I>N</I>, is twice the number of tissues.
 * <P>
 * The vector function is calculated as follows. For each <I>i</I>, 0 &lt;=
 * <I>i</I> &lt;= <I>M</I>&minus;1:
 * <CENTER>
 * <I>f</I><SUB><I>i</I></SUB>(<B>x</B>)&emsp;=&emsp;Model <I>S</I>(<I>t</I><SUB><I>i</I></SUB>;<B>x</B>) &minus; Measured <I>S</I>(<I>t</I><SUB><I>i</I></SUB>)
 * </CENTER>
 * <P>
 * The model spin signal is defined as follows:
 * <CENTER>
 * Model <I>S</I>(<I>t</I><SUB><I>i</I></SUB>;<B>x</B>)&emsp;=&emsp;&Sigma;<SUB><I>j</I></SUB>&emsp;<I>&rho;</I><SUB><I>j</I></SUB> [1 &minus; 2 exp(&minus;<I>R</I><SUB><I>j</I></SUB> <I>t</I><SUB><I>i</I></SUB>)]
 * </CENTER>
 * <P>
 * To find the model spin signal parameters <B>x</B> that best fit the measured
 * spin signal, use an instance of class SpinSignalDifference with the nonlinear
 * least squares algorithm in class {@linkplain
 * edu.rit.numeric.NonLinearLeastSquares NonLinearLeastSquares}.
 *
 * @author  Alan Kaminsky
 * @version 11-Jun-2008
 */
public class SpinSignalDifference
	implements VectorFunction
	{

// Hidden data members.

	private int M;
	private int L;
	private Series t_series;
	private Series S_series;

// Exported constructors.

	/**
	 * Construct a new spin signal difference function.
	 *
	 * @param  data_series  X-Y series of measured time values (X) and measured
	 *                      spin signal values (Y).
	 * @param  L            Number of tissues in the model spin signal.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>data_series</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the series' length is 0. Thrown if
	 *     <TT>L</TT> &lt;= 0.
	 */
	public SpinSignalDifference
		(XYSeries data_series,
		 int L)
		{
		this (data_series.xSeries(), data_series.ySeries(), L);
		}

	/**
	 * Construct a new spin signal difference function.
	 *
	 * @param  t_series  Series of measured time values.
	 * @param  S_series  Series of measured spin signal values.
	 * @param  L         Number of tissues in the model spin signal.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>t_series</TT> is null or
	 *     <TT>S_series</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if either series' length is 0. Thrown if
	 *     the two series have different lengths. Thrown if <TT>L</TT> &lt;= 0.
	 */
	public SpinSignalDifference
		(Series t_series,
		 Series S_series,
		 int L)
		{
		if (t_series == null)
			{
			throw new NullPointerException
				("SpinSignalDifference(): t_series is null");
			}
		if (S_series == null)
			{
			throw new NullPointerException
				("SpinSignalDifference(): S_series is null");
			}
		if (t_series.length() == 0)
			{
			throw new IllegalArgumentException
				("SpinSignalDifference(): t_series is zero length, illegal");
			}
		if (S_series.length() == 0)
			{
			throw new IllegalArgumentException
				("SpinSignalDifference(): S_series is zero length, illegal");
			}
		if (t_series.length() != S_series.length())
			{
			throw new IllegalArgumentException
				("SpinSignalDifference(): t_series and S_series have different lengths, illegal");
			}
		if (L <= 0)
			{
			throw new IllegalArgumentException
				("SpinSignalDifference(): L (= "+L+") <= 0, illegal");
			}

		this.M = t_series.length();
		this.L = L;
		this.t_series = t_series;
		this.S_series = S_series;
		}

// Exported operations.

	/**
	 * Returns the length of the result vector, <I>M</I>.
	 *
	 * @return  <I>M</I>.
	 */
	public int resultLength()
		{
		return M;
		}

	/**
	 * Returns the length of the argument vector, <I>N</I>.
	 *
	 * @return  <I>N</I>.
	 */
	public int argumentLength()
		{
		return L<<1;
		}

	/**
	 * Evaluate this function with the given argument vector. The result is
	 * stored in the vector <TT>y</TT>. Specifically, for <I>i</I> = 0 to
	 * <I>M</I>&minus;1, <I>y</I><SUB><I>i</I></SUB> =
	 * <I>f</I><SUB><I>i</I></SUB>(<B>x</B>).
	 *
	 * @param  x  Argument vector (input). Must be an <I>N</I>-element array.
	 * @param  y  Result vector (output). Must be an <I>M</I>-element array.
	 *
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if any argument in <TT>x</TT> is outside
	 *     the allowed set of values for this function.
	 * @exception  RangeException
	 *     (unchecked exception) Thrown if any element of the result vector is
	 *     outside the range of type <TT>double</TT>.
	 */
	public void f
		(double[] x,
		 double[] y)
		{
		// Subtract out the measured spin signal.
		for (int i = 0; i < M; ++ i)
			{
			y[i] = -S_series.x(i);
			}

		// Add in the model spin signal.
		for (int i = 0; i < M; ++ i)
			{
			double t_i = t_series.x(i);
			for (int j = 0; j < L; ++ j)
				{
				double dens_j = x[(j<<1)];
				double rate_j = x[(j<<1)+1];
				y[i] += SpinSignal.S (dens_j, rate_j, t_i);
				}
			}
		}

	/**
	 * Calculate this function's Jacobian matrix with the given argument vector.
	 * The result is stored in the matrix <TT>dydx</TT>. Specifically, for
	 * <I>i</I> = 0 to <I>M</I>&minus;1 and <I>j</I> = 0 to <I>N</I>&minus;1,
	 * <I>dydx</I><SUB><I>i,j</I></SUB> =
	 * &part;<I>f</I><SUB><I>i</I></SUB>(<B>x</B>)&nbsp;&frasl;&nbsp;&part;<I>x</I><SUB><I>j</I></SUB>&nbsp;.
	 *
	 * @param  x     Argument vector (input). Must be an <I>N</I>-element array.
	 * @param  dydx  Jacobian matrix (output). Must be an
	 *               <I>M</I>&times;<I>N</I>-element matrix.
	 *
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if any argument in <TT>x</TT> is outside
	 *     the allowed set of values for this function.
	 * @exception  RangeException
	 *     (unchecked exception) Thrown if any element of the Jacobian matrix is
	 *     outside the range of type <TT>double</TT>.
	 */
	public void df
		(double[] x,
		 double[][] dydx)
		{
		for (int i = 0; i < M; ++ i)
			{
			double t_i = t_series.x(i);
			double[] dydx_i = dydx[i];
			for (int j = 0; j < L; ++ j)
				{
				double dens_j = x[(j<<1)];
				double rate_j = x[(j<<1)+1];
				double twoexp = 2.0*Math.exp(-rate_j*t_i);
				dydx_i[(j<<1)] = 1.0 - twoexp;
				dydx_i[(j<<1)+1] = dens_j*t_i*twoexp;
				}
			}
		}

	}
