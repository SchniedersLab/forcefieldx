//******************************************************************************
//
// File:    Series.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.Series
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

import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Class Series is the abstract base class for a series of real values (type
 * <TT>double</TT>).
 *
 * @author  Alan Kaminsky
 * @version 25-Jul-2011
 */
public abstract class Series
	implements Iterable<Double>
	{

// Exported helper classes.

	/**
	 * Class Series.Stats holds the mean, variance, and standard deviation of a
	 * {@linkplain Series}.
	 * <P>
	 * If the series is empty, the mean, variance, and standard deviation are
	 * set to <TT>Double.NaN</TT> (not-a-number).
	 *
	 * @author  Alan Kaminsky
	 * @version 20-Jul-2011
	 */
	public static class Stats
		{
		/**
		 * Mean of the series' data values.
		 */
		public final double mean;

		/**
		 * Variance of the series' data values.
		 */
		public final double var;

		/**
		 * Standard deviation of the series' data values.
		 */
		public final double stddev;

		/**
		 * Construct a new Series.Stats object.
		 */
		private Stats
			(Series series)
			{
			int n = series.length();
			if (n == 0)
				{
				mean = Double.NaN;
				var = Double.NaN;
				stddev = Double.NaN;
				}
			else
				{
				double sum = 0.0;
				for (int i = 0; i < n; ++ i) sum += series.x(i);
				mean = sum/n;
				double sumdev = 0.0;
				double sumdevsqr = 0.0;
				for (int i = 0; i < n; ++ i)
					{
					double dev = series.x(i) - mean;
					sumdev += dev;
					sumdevsqr += dev*dev;
					}
				var = n == 1 ?  0.0 : (sumdevsqr - sumdev*sumdev/n)/(n - 1);
				stddev = Math.sqrt (var);
				}
			}
		}

	/**
	 * Class Series.RobustStats holds the median, mean absolute deviation, and
	 * quantiles of a {@linkplain Series}.
	 * <P>
	 * If the series is empty, the median, mean absolute deviation, and
	 * quantiles are set to <TT>Double.NaN</TT> (not-a-number), and the
	 * histogram bin counts are always 0.
	 *
	 * @author  Alan Kaminsky
	 * @version 25-Jul-2011
	 */
	public static class RobustStats
		{
		/**
		 * Median of the series' data values, equal to <TT>quantile(0.5)</TT>.
		 */
		public final double median;

		/**
		 * Mean absolute deviation of the series' data values. The absolute
		 * deviation of one value <I>x</I> is
		 * |<I>x</I>&nbsp;&minus;&nbsp;<I>median</I>|. The mean absolute
		 * deviation is the mean of all the absolute deviations.
		 */
		public final double meanAbsDev;

		/**
		 * The series' data values in ascending order.
		 */
		private double[] data;

		/**
		 * Construct a new Series.RobustStats object.
		 */
		private RobustStats
			(Series series)
			{
			// Sort series data.
			int n = series.length();
			data = new double [n];
			for (int i = 0; i < n; ++ i) data[i] = series.x(i);
			Arrays.sort (data);

			if (n == 0)
				{
				median = Double.NaN;
				meanAbsDev = Double.NaN;
				}
			else
				{
				// Compute median.
				median = quantile (0.5);

				// Compute mean absolute deviation.
				double sum = 0.0;
				for (int i = 0; i < n; ++ i)
					{
					sum += Math.abs (data[i] - median);
					}
				meanAbsDev = sum/n;
				}
			}

		/**
		 * Returns the given quantile of the series' data values. The
		 * <I>q</I>-th quantile is that value <I>x</I> such that a fraction
		 * <I>q</I> of the series' data values are less than or equal to
		 * <I>x</I>.
		 * <P>
		 * If the series is empty, then <TT>Double.NaN</TT> (not-a-number) is
		 * returned. If <I>q</I> is 0, then <TT>Double.NEGATIVE_INFINITY</TT> is
		 * returned.
		 *
		 * @param  q  Quantile, 0.0 &le; <I>q</I> &le; 1.0.
		 *
		 * @return  The <I>q</I>-th quantile of the series' data values.
		 *
		 * @exception  IllegalArgumentException
		 *     (unchecked exception) Thrown if <I>q</I> is out of range.
		 */
		public double quantile
			(double q)
			{
			if (0.0 > q || q > 1.0)
				{
				throw new IllegalArgumentException
					("Series.RobustStats.quantile(): q = "+q+" illegal");
				}
			return
				data.length == 0 ?
					Double.NaN :
				q == 0.0 ?
					Double.NEGATIVE_INFINITY :
					data[(int)(Math.ceil (q*data.length)) - 1];
			}

		/**
		 * Returns the count for the given histogram bin. This is the number of
		 * data values in the series that are greater than or equal to
		 * <TT>lb</TT> and less than <TT>ub</TT>.
		 *
		 * @param  lb  Lower bound (inclusive) of histogram bin.
		 * @param  ub  Upper bound (exclusive) of histogram bin.
		 *
		 * @return  Histogram bin count.
		 *
		 * @exception  IllegalArgumentException
		 *     (unchecked exception) Thrown if <TT>lb</TT> &gt; <TT>ub</TT>.
		 */
		public int histogram
			(double lb,
			 double ub)
			{
			if (lb > ub)
				{
				throw new IllegalArgumentException
					("RobustStats.histogram(): lb ("+lb+") > ub ("+ub+
					 ") illegal");
				}
			int n = data.length;
			int i = 0;
			while (i < n && data[i] < lb) ++ i;
			int j = i;
			while (j < n && data[j] < ub) ++ j;
			return j - i;
			}
		}

// Exported constructors.

	/**
	 * Construct a new series.
	 */
	public Series()
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
	 * Returns the minimum value in this series.
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
	 * Returns the maximum value in this series.
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
	 * Returns a {@linkplain Stats Stats} object containing statistics of this
	 * series.
	 * <P>
	 * <I>Note:</I> The returned object contains the statistics of a
	 * <I>snapshot</I> of this series at the time <TT>stats()</TT> was called.
	 * Changing the data in this series will <I>not</I> change the contents of
	 * the returned object.
	 *
	 * @return  Statistics.
	 */
	public Stats stats()
		{
		return new Stats (this);
		}

	/**
	 * Returns a {@linkplain RobustStats RobustStats} object containing robust
	 * statistics of this series.
	 * <P>
	 * <I>Note:</I> The returned object contains the statistics of a
	 * <I>snapshot</I> of this series at the time <TT>robustStats()</TT> was
	 * called. Changing the data in this series will <I>not</I> change the
	 * contents of the returned object.
	 *
	 * @return  Robust statistics.
	 */
	public RobustStats robustStats()
		{
		return new RobustStats (this);
		}

	/**
	 * Returns an iterator over the values in this series.
	 *
	 * @return  Iterator.
	 */
	public Iterator<Double> iterator()
		{
		return new Iterator<Double>()
			{
			int i = 0;

			public boolean hasNext()
				{
				return i < length();
				}

			public Double next()
				{
				try
					{
					return x (i ++);
					}
				catch (ArrayIndexOutOfBoundsException exc)
					{
					throw new NoSuchElementException();
					}
				}

			public void remove()
				{
				throw new UnsupportedOperationException();
				}
			};
		}

	/**
	 * Print this series on the standard output. Each line of output consists of
	 * the index and the value, separated by a tab.
	 */
	public void print()
		{
		print (System.out);
		}

	/**
	 * Print this series on the given print stream. Each line of output consists
	 * of the index and the value, separated by a tab.
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
			theStream.println (x(i));
			}
		}

	/**
	 * Print this series on the given print writer. Each line of output consists
	 * of the index and the value, separated by a tab.
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
			theWriter.println (x(i));
			}
		}

	}
