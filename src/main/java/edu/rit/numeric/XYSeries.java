//******************************************************************************
//
// File:    XYSeries.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.XYSeries
//
// This Java source file is copyright (C) 2011 by Alan Kaminsky. All rights
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

import java.io.PrintStream;
import java.io.PrintWriter;

/**
 * Class XYSeries is the abstract base class for a series of (<I>x,y</I>) pairs
 * of real values (type <TT>double</TT>).
 *
 * @author  Alan Kaminsky
 * @version 20-Jul-2011
 */
public abstract class XYSeries
	{

// Exported helper classes.

	/**
	 * Class XYSeries.Regression holds the results of a regression on an
	 * {@linkplain XYSeries}.
	 *
	 * @author  Alan Kaminsky
	 * @version 12-Jun-2007
	 */
	public static class Regression
		{
		/**
		 * Intercept <I>a</I>.
		 */
		public final double a;

		/**
		 * Slope <I>b</I>.
		 */
		public final double b;

		/**
		 * Correlation.
		 */
		public final double corr;

		/**
		 * Construct a new Regression object.
		 */
		private Regression
			(double a,
			 double b,
			 double corr)
			{
			this.a = a;
			this.b = b;
			this.corr = corr;
			}
		}

	/**
	 * Class XYSeries.XSeriesView provides a series view of the X values in an
	 * XY series.
	 *
	 * @author  Alan Kaminsky
	 * @version 12-Jun-2007
	 */
	private static class XSeriesView
		extends Series
		{
		private XYSeries outer;

		public XSeriesView
			(XYSeries outer)
			{
			this.outer = outer;
			}

		public int length()
			{
			return outer.length();
			}

		public double x
			(int i)
			{
			return outer.x (i);
			}
		}

	/**
	 * Class XYSeries.YSeriesView provides a series view of the Y values in an
	 * XY series.
	 *
	 * @author  Alan Kaminsky
	 * @version 12-Jun-2007
	 */
	private static class YSeriesView
		extends Series
		{
		private XYSeries outer;

		public YSeriesView
			(XYSeries outer)
			{
			this.outer = outer;
			}

		public int length()
			{
			return outer.length();
			}

		public double x
			(int i)
			{
			return outer.y (i);
			}
		}

// Exported constructors.

	/**
	 * Construct a new XY series.
	 */
	public XYSeries()
		{
		}

// Exported operations.

	/**
	 * Returns the number of values in this series.
	 *
	 * @return  Length.
	 */
	public abstract int length();

	/**
	 * Determine if this series is empty.
	 *
	 * @return  True if this series is empty (length = 0), false otherwise.
	 */
	public boolean isEmpty()
		{
		return length() == 0;
		}

	/**
	 * Returns the given X value in this series.
	 *
	 * @param  i  Index.
	 *
	 * @return  The X value in this series at index <TT>i</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range
	 *     <TT>0</TT> .. <TT>length()-1</TT>.
	 */
	public abstract double x
		(int i);

	/**
	 * Returns the given Y value in this series.
	 *
	 * @param  i  Index.
	 *
	 * @return  The Y value in this series at index <TT>i</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range
	 *     <TT>0</TT> .. <TT>length()-1</TT>.
	 */
	public abstract double y
		(int i);

	/**
	 * Returns the minimum X value in this series.
	 *
	 * @return  Minimum X value.
	 */
	public double minX()
		{
		int n = length();
		double result = Double.POSITIVE_INFINITY;
		for (int i = 0; i < n; ++ i) result = Math.min (result, x(i));
		return result;
		}

	/**
	 * Returns the maximum X value in this series.
	 *
	 * @return  Maximum X value.
	 */
	public double maxX()
		{
		int n = length();
		double result = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < n; ++ i) result = Math.max (result, x(i));
		return result;
		}

	/**
	 * Returns the minimum Y value in this series.
	 *
	 * @return  Minimum Y value.
	 */
	public double minY()
		{
		int n = length();
		double result = Double.POSITIVE_INFINITY;
		for (int i = 0; i < n; ++ i) result = Math.min (result, y(i));
		return result;
		}

	/**
	 * Returns the maximum Y value in this series.
	 *
	 * @return  Maximum Y value.
	 */
	public double maxY()
		{
		int n = length();
		double result = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < n; ++ i) result = Math.max (result, y(i));
		return result;
		}

	/**
	 * Returns the linear regression of the (<I>x,y</I>) values in this XY
	 * series. The linear function <I>y</I> = <I>a</I> + <I>bx</I> is fitted to
	 * the data. The return value is a {@linkplain Regression Regression} object
	 * containing the intercept <I>a,</I> the slope <I>b,</I> and the
	 * correlation, respectively.
	 * <P>
	 * <I>Note:</I> The returned object contains the regression of a
	 * <I>snapshot</I> of this series at the time <TT>linearRegression()</TT>
	 * was called. Changing the data in this series will <I>not</I> change the
	 * contents of the returned object.
	 *
	 * @return  Regression.
	 */
	public Regression linearRegression()
		{
		// Accumulate sums.
		int n = length();
		double sum_x = 0.0;
		double sum_y = 0.0;
		for (int i = 0; i < n; ++ i)
			{
			sum_x += x(i);
			sum_y += y(i);
			}

		// Compute means of X and Y.
		double mean_x = sum_x / n;
		double mean_y = sum_y / n;

		// Compute variances of X and Y.
		double xt;
		double yt;
		double sum_xt_sqr = 0.0;
		double sum_yt_sqr = 0.0;
		double sum_xt_yt  = 0.0;
		double b = 0.0;
		for (int i = 0; i < n; ++ i)
			{
			xt = x(i) - mean_x;
			yt = y(i) - mean_y;
			sum_xt_sqr += xt * xt;
			sum_yt_sqr += yt * yt;
			sum_xt_yt  += xt * yt;
			b += xt * y(i);
			}

		// Compute results.
		b /= sum_xt_sqr;
		double a = (sum_y - sum_x * b) / n;
		double corr = sum_xt_yt / (Math.sqrt (sum_xt_sqr * sum_yt_sqr) + TINY);
		return new Regression (a, b, corr);
		}

	private static final double TINY = 1.0e-20;

	/**
	 * Returns a {@linkplain Series} view of the X values in this XY series.
	 * <P>
	 * <I>Note:</I> The returned Series object is backed by this XY series
	 * object. Changing the contents of this XY series object will change the
	 * contents of the returned Series object.
	 *
	 * @return  Series of X values.
	 */
	public Series xSeries()
		{
		return new XSeriesView (this);
		}

	/**
	 * Returns a {@linkplain Series} view of the Y values in this XY series.
	 * <P>
	 * <I>Note:</I> The returned Series object is backed by this XY series
	 * object. Changing the contents of this XY series object will change the
	 * contents of the returned Series object.
	 *
	 * @return  Series of Y values.
	 */
	public Series ySeries()
		{
		return new YSeriesView (this);
		}

	/**
	 * Print this XY series on the standard output. Each line of output consists
	 * of the index, the <I>x</I> value, and the <I>y</I> value, separated by
	 * tabs.
	 */
	public void print()
		{
		print (System.out);
		}

	/**
	 * Print this XY series on the given print stream. Each line of output
	 * consists of the index, the <I>x</I> value, and the <I>y</I> value,
	 * separated by tabs.
	 *
	 * @param  theStream  Print stream.
	 */
	public void print
		(PrintStream theStream)
		{
		int n = length();
		for (int i = 0; i < n; ++ i)
			{
			theStream.print (i);
			theStream.print ('\t');
			theStream.print (x(i));
			theStream.print ('\t');
			theStream.println (y(i));
			}
		}

	/**
	 * Print this XY series on the given print writer. Each line of output
	 * consists of the index, the <I>x</I> value, and the <I>y</I> value,
	 * separated by tabs.
	 *
	 * @param  theWriter  Print writer.
	 */
	public void print
		(PrintWriter theWriter)
		{
		int n = length();
		for (int i = 0; i < n; ++ i)
			{
			theWriter.print (i);
			theWriter.print ('\t');
			theWriter.print (x(i));
			theWriter.print ('\t');
			theWriter.println (y(i));
			}
		}

	}
