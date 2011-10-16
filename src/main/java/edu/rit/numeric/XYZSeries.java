//******************************************************************************
//
// File:    XYZSeries.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.XYZSeries
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
 * Class XYZSeries is the abstract base class for a series of (<I>x,y,z</I>)
 * triples of real values (type <TT>double</TT>).
 *
 * @author  Alan Kaminsky
 * @version 20-Jul-2011
 */
public abstract class XYZSeries
	{

// Exported helper classes.

	/**
	 * Class XYZSeries.Regression holds the results of a regression on an
	 * {@linkplain XYZSeries}. The <I>x</I> and <I>y</I> values are the data
	 * points; each <I>z</I> value is the standard deviation of the measurement
	 * error in the corresponding <I>y</I> value, assuming the measurement error
	 * distribution is a zero-mean Gaussian.
	 *
	 * @author  Alan Kaminsky
	 * @version 16-Mar-2010
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
		 * Variance in the estimate of <I>a</I>.
		 */
		public final double var_a;

		/**
		 * Variance in the estimate of <I>b</I>.
		 */
		public final double var_b;

		/**
		 * Covariance in the estimates of <I>a</I> and <I>b</I>.
		 */
		public final double cov_ab;

		/**
		 * Coefficient of correlation between the uncertainty in <I>a</I> and
		 * the uncertainty in <I>b</I>.
		 */
		public final double r_ab;

		/**
		 * Chi-square merit function of the fit between the model (this
		 * Regression object) and the data (original XYZSeries object).
		 */
		public final double chi2;

		/**
		 * Significance of <I>&chi;</I><SUP>2</SUP>.
		 */
		public final double significance;

		/**
		 * Construct a new Regression object.
		 */
		private Regression
			(double a,
			 double b,
			 double var_a,
			 double var_b,
			 double cov_ab,
			 double r_ab,
			 double chi2,
			 double significance)
			{
			this.a = a;
			this.b = b;
			this.var_a = var_a;
			this.var_b = var_b;
			this.cov_ab = cov_ab;
			this.r_ab = r_ab;
			this.chi2 = chi2;
			this.significance = significance;
			}
		}

	/**
	 * Class XYZSeries.XSeriesView provides a series view of the X values in an
	 * XYZ series.
	 *
	 * @author  Alan Kaminsky
	 * @version 12-Jun-2007
	 */
	private static class XSeriesView
		extends Series
		{
		private XYZSeries outer;

		public XSeriesView
			(XYZSeries outer)
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
	 * Class XYZSeries.YSeriesView provides a series view of the Y values in an
	 * XYZ series.
	 *
	 * @author  Alan Kaminsky
	 * @version 12-Jun-2007
	 */
	private static class YSeriesView
		extends Series
		{
		private XYZSeries outer;

		public YSeriesView
			(XYZSeries outer)
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

	/**
	 * Class XYZSeries.ZSeriesView provides a series view of the Z values in an
	 * XYZ series.
	 *
	 * @author  Alan Kaminsky
	 * @version 12-Jun-2007
	 */
	private static class ZSeriesView
		extends Series
		{
		private XYZSeries outer;

		public ZSeriesView
			(XYZSeries outer)
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
			return outer.z (i);
			}
		}

	/**
	 * Class XYZSeries.XYSeriesView provides an XY series view of the X and Y
	 * values in an XYZ series.
	 *
	 * @author  Alan Kaminsky
	 * @version 13-Oct-2007
	 */
	private static class XYSeriesView
		extends XYSeries
		{
		private XYZSeries outer;

		public XYSeriesView
			(XYZSeries outer)
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

		public double y
			(int i)
			{
			return outer.y (i);
			}
		}

// Exported constructors.

	/**
	 * Construct a new XYZ series.
	 */
	public XYZSeries()
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
	 * Returns the given Z value in this series.
	 *
	 * @param  i  Index.
	 *
	 * @return  The Z value in this series at index <TT>i</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range
	 *     <TT>0</TT> .. <TT>length()-1</TT>.
	 */
	public abstract double z
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
		for (int i = 0; i < n; ++ i)
			{
			result = Math.min (result, x(i));
			}
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
		for (int i = 0; i < n; ++ i)
			{
			result = Math.max (result, x(i));
			}
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
		for (int i = 0; i < n; ++ i)
			{
			result = Math.min (result, y(i));
			}
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
		for (int i = 0; i < n; ++ i)
			{
			result = Math.max (result, y(i));
			}
		return result;
		}

	/**
	 * Returns the minimum Z value in this series.
	 *
	 * @return  Minimum Z value.
	 */
	public double minZ()
		{
		int n = length();
		double result = Double.POSITIVE_INFINITY;
		for (int i = 0; i < n; ++ i)
			{
			result = Math.min (result, z(i));
			}
		return result;
		}

	/**
	 * Returns the maximum Z value in this series.
	 *
	 * @return  Maximum Z value.
	 */
	public double maxZ()
		{
		int n = length();
		double result = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < n; ++ i)
			{
			result = Math.max (result, z(i));
			}
		return result;
		}

	/**
	 * Returns the linear regression of the (<I>x,y</I>) values in this XYZ
	 * series. Each <I>z</I> value is the standard deviation of the measurement
	 * error in the corresponding <I>y</I> value, assuming the measurement error
	 * distribution is a zero-mean Gaussian. The linear function
	 * <I>y</I>&nbsp;=&nbsp;<I>a</I>&nbsp;+&nbsp;<I>bx</I> is fitted to the
	 * data. The return value is a {@linkplain Regression Regression} object
	 * containing the intercept <I>a,</I> the slope <I>b,</I> and various
	 * statistics.
	 *
	 * @return  Regression.
	 */
	public Regression linearRegression()
		{
		// Accumulate sums.
		int n = length();
		double Sx = 0.0;
		double Sy = 0.0;
		double S = 0.0;
		for (int i = 0; i < n; ++ i)
			{
			double sigma_i = z(i);
			double sigma2_i = sigma_i*sigma_i;
			Sx += x(i)/sigma2_i;
			Sy += y(i)/sigma2_i;
			S += 1.0/sigma2_i;
			}
		double Sx_over_S = Sx/S;

		// Calculate b.
		double b = 0.0;
		double Stt = 0.0;
		for (int i = 0; i < n; ++ i)
			{
			double sigma_i = z(i);
			double t_i = (x(i) - Sx_over_S)/sigma_i;
			b += t_i*y(i)/sigma_i;
			Stt += t_i*t_i;
			}
		b /= Stt;

		// Calculate other model parameters.
		double a = (Sy - Sx*b)/S;
		double var_a = (1.0 + Sx*Sx/S/Stt)/S;
		double var_b = 1.0/Stt;
		double cov_ab = -Sx/S/Stt;
		double r_ab = cov_ab/Math.sqrt(var_a)/Math.sqrt(var_b);

		// Calculate chi^2.
		double chi2 = 0.0;
		for (int i = 0; i < n; ++ i)
			{
			double d = (y(i) - a - b*x(i))/z(i);
			chi2 += d*d;
			}
		double significance = Statistics.chiSquarePvalue (n - 2, chi2);

		// Return results.
		return new Regression
			(a, b, var_a, var_b, cov_ab, r_ab, chi2, significance);
		}

	/**
	 * Returns a {@linkplain Series} view of the X values in this XYZ series.
	 * <P>
	 * <I>Note:</I> The returned Series object is backed by this XYZ series
	 * object. Changing the contents of this XYZ series object will change the
	 * contents of the returned Series object.
	 *
	 * @return  Series of X values.
	 */
	public Series xSeries()
		{
		return new XSeriesView (this);
		}

	/**
	 * Returns a {@linkplain Series} view of the Y values in this XYZ series.
	 * <P>
	 * <I>Note:</I> The returned Series object is backed by this XYZ series
	 * object. Changing the contents of this XYZ series object will change the
	 * contents of the returned Series object.
	 *
	 * @return  Series of Y values.
	 */
	public Series ySeries()
		{
		return new YSeriesView (this);
		}

	/**
	 * Returns a {@linkplain Series} view of the Z values in this XYZ series.
	 * <P>
	 * <I>Note:</I> The returned Series object is backed by this XYZ series
	 * object. Changing the contents of this XYZ series object will change the
	 * contents of the returned Series object.
	 *
	 * @return  Series of Z values.
	 */
	public Series zSeries()
		{
		return new ZSeriesView (this);
		}

	/**
	 * Returns an {@linkplain XYSeries} view of the X and Y values in this XYZ
	 * series.
	 * <P>
	 * <I>Note:</I> The returned XYSeries object is backed by this XYZ series
	 * object. Changing the contents of this XYZ series object will change the
	 * contents of the returned Series object.
	 *
	 * @return  Series of X and Y values.
	 */
	public XYSeries xySeries()
		{
		return new XYSeriesView (this);
		}

	/**
	 * Print this XYZ series on the standard output. Each line of output
	 * consists of the index, the <I>x</I> value, the <I>y</I> value, and the
	 * <I>z</I> value, separated by tabs.
	 */
	public void print()
		{
		print (System.out);
		}

	/**
	 * Print this XYZ series on the given print stream. Each line of output
	 * consists of the index, the <I>x</I> value, the <I>y</I> value, and the
	 * <I>z</I> value, separated by tabs.
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
			theStream.print (y(i));
			theStream.print ('\t');
			theStream.println (z(i));
			}
		}

	/**
	 * Print this XYZ series on the given print writer. Each line of output
	 * consists of the index, the <I>x</I> value, the <I>y</I> value, and the
	 * <I>z</I> value, separated by tabs.
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
			theWriter.print (y(i));
			theWriter.print ('\t');
			theWriter.println (z(i));
			}
		}

	}
