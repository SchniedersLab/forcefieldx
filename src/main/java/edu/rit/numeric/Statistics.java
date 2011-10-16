//******************************************************************************
//
// File:    Statistics.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.Statistics
//
// This Java source file is copyright (C) 2009 by Alan Kaminsky. All rights
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

import java.math.BigInteger;

// For unit test main program.
//import edu.rit.util.Random;
//import java.util.Arrays;

/**
 * Class Statistics provides static methods for doing statistical tests.
 * <P>
 * For each statistical test, there is a method that returns the "p-value" of
 * the test statistic. This is the probability that the test statistic would
 * have a value greater than or equal to the observed value if the null
 * hypothesis is true.
 *
 * @author  Alan Kaminsky
 * @version 28-Apr-2010
 */
public class Statistics
	{

// Prevent construction.

	private Statistics()
		{
		}

// Exported operations.

	/**
	 * Do a chi-square test on the given data. The null hypothesis is that the
	 * data was drawn from the distribution given by <TT>expected</TT>. The
	 * <TT>measured</TT> and <TT>expected</TT> arrays must be the same length.
	 *
	 * @param  measured  Measured count in each bin.
	 * @param  expected  Expected count in each bin.
	 *
	 * @return  Chi-square statistic.
	 */
	public static double chiSquareTest
		(double[] measured,
		 double[] expected)
		{
		double chisqr = 0.0;
		for (int i = 0; i < measured.length; ++ i)
			{
			double d = measured[i] - expected[i];
			chisqr += d*d/expected[i];
			}
		return chisqr;
		}

	/**
	 * Returns the p-value of a chi-square statistic.
	 *
	 * @param  N       Degrees of freedom.
	 * @param  chisqr  Chi-square statistic.
	 *
	 * @return  P-value.
	 */
	public static double chiSquarePvalue
		(double N,
		 double chisqr)
		{
		return gammq (0.5*N, 0.5*chisqr);
		}

	/**
	 * Do a Bernoulli chi-square test on the given data. The null hypothesis is
	 * that the data was drawn from a Bernoulli distribution with both outcomes
	 * equally likely (e.g., a fair coin). <TT>total</TT> is the total number of
	 * trials. <TT>measured</TT> is the number of trials yielding one of the
	 * outcomes. (<TT>total</TT>&nbsp;&minus;&nbsp;<TT>measured</TT>) is the
	 * number of trials yielding the other outcome.
	 *
	 * @param  total     Total number of trials.
	 * @param  measured  Number of trials yielding one of the outcomes.
	 *
	 * @return  Chi-square statistic.
	 */
	public static double bernoulliChiSquareTest
		(long total,
		 long measured)
		{
		double expected = 0.5*total;
		double d = measured - expected;
		return 2.0*d*d/expected;
		}

	/**
	 * Returns the p-value of a Bernoulli chi-square statistic.
	 *
	 * @param  chisqr  Chi-square statistic.
	 *
	 * @return  P-value.
	 */
	public static double bernoulliChiSquarePvalue
		(double chisqr)
		{
		return gammq (0.5, 0.5*chisqr);
		}

	/**
	 * Do a Y-square test on the given data. The null hypothesis is that the
	 * data was drawn from the distribution given by <TT>expected</TT>. The
	 * <TT>measured</TT> and <TT>expected</TT> arrays must be the same length.
	 * <P>
	 * The Y-square test is similar to the chi-square test, except the Y-square
	 * statistic is valid even if the expected counts in some of the bins are
	 * small, which is not true of the chi-square statistic. For further
	 * information, see:
	 * <P>
	 * L. Lucy. Hypothesis testing for meagre data sets. <I>Monthly Notices of
	 * the Royal Astronomical Society,</I> 318(1):92-100, October 2000.
	 *
	 * @param  N         Degrees of freedom.
	 * @param  measured  Measured count in each bin.
	 * @param  expected  Expected count in each bin.
	 *
	 * @return  Y-square statistic.
	 */
	public static double ySquareTest
		(int N,
		 double[] measured,
		 double[] expected)
		{
		double twoN = 2.0*N;
		double sum = 0.0;
		for (int i = 0; i < expected.length; ++ i)
			{
			sum += 1.0/expected[i];
			}
		return N + Math.sqrt(twoN/(twoN + sum))*
					(chiSquareTest (measured, expected) - N);
		}

	/**
	 * Returns the p-value of a Y-square statistic.
	 *
	 * @param  N     Degrees of freedom.
	 * @param  ysqr  Y-square statistic.
	 *
	 * @return  P-value.
	 */
	public static double ySquarePvalue
		(double N,
		 double ysqr)
		{
		return gammq (0.5*N, 0.5*ysqr);
		}

	/**
	 * Do a Kolmogorov-Smirnov (K-S) test on the given data. The null hypothesis
	 * is that the data was drawn from a uniform distribution between 0.0 and
	 * 1.0.
	 * <P>
	 * The values in the <TT>data</TT> array must all be in the range 0.0
	 * through 1.0 and must be in ascending numerical order. The
	 * <TT>ksTest()</TT> method does not sort the data itself because the
	 * process that produced the data might already have sorted the data. If
	 * necessary, call <TT>Arrays.sort(data)</TT> before calling
	 * <TT>ksTest(data)</TT>.
	 *
	 * @param  data  Data array.
	 *
	 * @return  K-S statistic.
	 */
	public static double ksTest
		(double[] data)
		{
		int M = data.length;
		double N = M;
		double D = 0.0;
		double F_lower = 0.0;
		double F_upper;
		for (int i = 0; i < M; ++ i)
			{
			F_upper = (i+1) / N;
			D = Math.max (D, Math.abs (data[i] - F_lower));
			D = Math.max (D, Math.abs (data[i] - F_upper));
			F_lower = F_upper;
			}
		return D;
		}

	/**
	 * Do a Kolmogorov-Smirnov (K-S) test on the given data. The null hypothesis
	 * is that the data was drawn from the distribution specified by the given
	 * {@linkplain Function}. <TT>cdf.f(x)</TT> must return the value of the
	 * cumulative distribution function at <I>x</I>, in the range 0.0 through
	 * 1.0.
	 * <P>
	 * The values in the <TT>data</TT> array must all be in the domain of
	 * <TT>cdf</TT> and must be in ascending numerical order. The
	 * <TT>ksTest()</TT> method does not sort the data itself because the
	 * process that produced the data might already have sorted the data. If
	 * necessary, call <TT>Arrays.sort(data)</TT> before calling
	 * <TT>ksTest(data,cdf)</TT>.
	 *
	 * @param  data  Data array.
	 * @param  cdf   Cumulative distribution function.
	 *
	 * @return  K-S statistic.
	 */
	public static double ksTest
		(double[] data,
		 Function cdf)
		{
		int M = data.length;
		double N = M;
		double D = 0.0;
		double F_lower = 0.0;
		double F_upper;
		double cdf_i;
		for (int i = 0; i < M; ++ i)
			{
			F_upper = (i+1) / N;
			cdf_i = cdf.f (data[i]);
			D = Math.max (D, Math.abs (cdf_i - F_lower));
			D = Math.max (D, Math.abs (cdf_i - F_upper));
			F_lower = F_upper;
			}
		return D;
		}

	/**
	 * Do a Kolmogorov-Smirnov (K-S) test on the given data. The null hypothesis
	 * is that the data was drawn from a binomial random variable <I>X</I> that
	 * is the sum of <I>n</I> equiprobable Bernoulli random variables. For 0
	 * &le; <I>k</I> &le; <I>n</I>, the probability that <I>X</I> equals
	 * <I>k</I> is
	 * <P>
	 * <CENTER>
	 * Pr[<I>X</I> = <I>k</I>] = 2<SUP>&minus;<I>n</I></SUP> <I>n</I>! / <I>k</I>! / (<I>n</I> &minus; <I>k</I>)!
	 * </CENTER>
	 * <P>
	 * The values in the <TT>data</TT> array must all be in the range 0 ..
	 * <I>n</I> and must be in ascending numerical order. The
	 * <TT>binomialKsTest()</TT> method does not sort the data itself because
	 * the process that produced the data might already have sorted the data. If
	 * necessary, call <TT>Arrays.sort(data)</TT> before calling
	 * <TT>binomialKsTest(data,n)</TT>.
	 * <P>
	 * <I>Note:</I> To prevent roundoff error, the internal calculations are
	 * done using exact rational arithmetic. The final K-S statistic is then
	 * converted to a double-precision floating-point number and is returned.
	 *
	 * @param  data  Data array.
	 * @param  n     Number of Bernoulli random variables.
	 *
	 * @return  K-S statistic.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>n</TT> &le; 0.
	 */
	public static double binomialKsTest
		(int[] data,
		 int n)
		{
		if (n <= 0L)
			{
			throw new IllegalArgumentException
				("Statistics.binomialKsTest(): n = "+n+" illegal");
			}
		int M = data.length;
		BigRational D = new BigRational (0);
		BigRational F_lower = new BigRational();
		BigRational F_upper = new BigRational();
		int pf_i = 0;
		BigInteger pf_numer = new BigInteger ("1");
		BigInteger pf_numer_sum = new BigInteger ("0");
		BigInteger pf_denom = new BigInteger ("2") .pow (n);
		BigRational cdf = new BigRational (0);
		int i = 0;
		int j = 0;
		while (i < M)
			{
			j = i;
			while (j+1 < M && data[j+1] == data[i]) ++ j;
			F_lower.assign (i, M);
			F_upper.assign (j+1, M);
			while (pf_i < data[i])
				{
				pf_numer = pf_numer.multiply (BigInteger.valueOf (n - pf_i));
				++ pf_i;
				pf_numer = pf_numer.divide (BigInteger.valueOf (pf_i));
				pf_numer_sum = pf_numer_sum.add (pf_numer);
				}
			cdf.assign (pf_numer_sum, pf_denom);
			F_lower.sub (cdf) .abs();
			F_upper.sub (cdf) .abs();
			D.max (F_lower) .max (F_upper) .normalize();
			i = j+1;
			}
		return D.doubleValue();
		}

	/**
	 * Returns the p-value of a K-S statistic.
	 *
	 * @param  N  Number of data points.
	 * @param  D  K-S statistic.
	 *
	 * @return  P-value.
	 */
	public static double ksPvalue
		(double N,
		 double D)
		{
		double sqrt_N = Math.sqrt(N);
		double x = (sqrt_N + 0.12 + 0.11/sqrt_N) * D;
		x = -2.0*x*x;
		double a = 2.0;
		double sum = 0.0;
		double term;
		double absterm;
		double prevterm = 0.0;
		for (int j = 1; j <= 100; ++ j)
			{
			term = a * Math.exp (x*j*j);
			sum += term;
			absterm = Math.abs(term);
			if (absterm <= 1.0e-6*prevterm || absterm <= 1.0e-12*sum)
				{
				return sum;
				}
			a = -a;
			prevterm = absterm;
			}
		return 1.0; // Failed to converge
		}

	/**
	 * Returns the p-value of a statistic drawn from a normal distribution.
	 *
	 * @param  x       Statistic.
	 * @param  mean    Mean of the normal distribution.
	 * @param  stddev  Standard deviation of the normal distribution.
	 *
	 * @return  P-value.
	 */
	public static double normalPvalue
		(double x,
		 double mean,
		 double stddev)
		{
		return 1.0 - erfc (-INV_SQRT_2*(x - mean)/stddev);
		}

// Hidden operations.

	private static final double INV_SQRT_2 = 1.0/Math.sqrt(2.0);

	private static final double[] LGAMMA_COF = new double[]
		{ 57.1562356658629235,
		 -59.5979603554754912,
		  14.1360979747417471,
		  -0.491913816097620199,
		   0.339946499848118887e-4,
		   0.465236289270485756e-4,
		  -0.983744753048795646e-4,
		   0.158088703224912494e-3,
		  -0.210264441724104883e-3,
		   0.217439618115212643e-3,
		  -0.164318106536763890e-3,
		   0.844182239838527433e-4,
		  -0.261908384015814087e-4,
		   0.368991826595316234e-5};

	/**
	 * Returns the log-gamma function ln(gamma(x)). Assumes <TT>x</TT> &gt; 0.
	 */
	private static double lgamma
		(double x)
		{
		double y, tmp, ser;
		int i;

		y = x;
		tmp = x + 5.2421875;
		tmp = (x + 0.5)*Math.log(tmp) - tmp;
		ser = 0.999999999999997092;
		for (i = 0; i < LGAMMA_COF.length; ++ i)
			{
			y += 1.0;
			ser += LGAMMA_COF[i]/y;
			}
		return tmp + Math.log(2.5066282746310005*ser/x);
		}

	private static final int GAMMA_ITMAX = 200;
	private static final double GAMMA_EPS = 2.22e-16;
	private static final double GAMMA_FPMIN = (2.23e-308/GAMMA_EPS);

	/**
	 * Returns the incomplete gamma function P(a,x), evaluated by its series
	 * representation. Assumes <TT>a</TT> &gt; 0 and <TT>x</TT> &ge; 0.
	 */
	private static double gser
		(double a,
		 double x)
		{
		double ap, del, sum;
		int i;

		ap = a;
		del = 1.0/a;
		sum = del;
		for (i = 1; i <= GAMMA_ITMAX; ++ i)
			{
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if (Math.abs(del) < Math.abs(sum)*GAMMA_EPS)
				{
				return sum*Math.exp(-x + a*Math.log(x) - lgamma(a));
				}
			}
		return 1.0; // Too many iterations
		}

	/**
	 * Returns the complementary incomplete gamma function Q(a,x), evaluated by
	 * its continued fraction representation. Assumes <TT>a</TT> &gt; 0 and
	 * <TT>x</TT> &ge; 0.
	 */
	private static double gcf
		(double a,
		 double x)
		{
		double b, c, d, h, an, del;
		int i;

		b = x + 1.0 - a;
		c = 1.0/GAMMA_FPMIN;
		d = 1.0/b;
		h = d;
		for (i = 1; i <= GAMMA_ITMAX; ++ i)
			{
			an = -i*(i - a);
			b += 2.0;
			d = an*d + b;
			if (Math.abs(d) < GAMMA_FPMIN) d = GAMMA_FPMIN;
			c = b + an/c;
			if (Math.abs(c) < GAMMA_FPMIN) c = GAMMA_FPMIN;
			d = 1.0/d;
			del = d*c;
			h *= del;
			if (Math.abs(del - 1.0) < GAMMA_EPS)
				{
				return Math.exp(-x + a*Math.log(x) - lgamma(a))*h;
				}
			}
		return 0.0; // Too many iterations
		}

	/**
	 * Returns the incomplete gamma function P(a,x).
	 */
	private static double gammp
		(double a,
		 double x)
		{
		if (a <= 0.0)
			{
			throw new IllegalArgumentException ("gammp(): a = "+a+" illegal");
			}
		if (x < 0.0)
			{
			throw new IllegalArgumentException ("gammp(): x = "+x+" illegal");
			}
		return x == 0.0 ? 0.0 : x < a + 1.0 ? gser(a,x) : 1.0 - gcf(a,x);
		}

	/**
	 * Returns the complementary incomplete gamma function Q(a,x) = 1 - P(a,x).
	 */
	private static double gammq
		(double a,
		 double x)
		{
		if (a <= 0.0)
			{
			throw new IllegalArgumentException ("gammq(): a = "+a+" illegal");
			}
		if (x < 0.0)
			{
			throw new IllegalArgumentException ("gammq(): x = "+x+" illegal");
			}
		return x == 0.0 ? 1.0 : x < a + 1.0 ? 1.0 - gser(a,x) : gcf(a,x);
		}

	/**
	 * Returns the complementary error function erfc(x).
	 */
	private static double erfc
		(double x)
		{
		return gammq (0.5, x*x);
		}

// Unit test main program.

//	/**
//	 * Unit test main program. Does a K-S test on N random doubles, prints the
//	 * K-S statistic, and prints the p-value.
//	 * <P>
//	 * Usage: java edu.rit.numeric.Statistics <I>seed</I> <I>N</I>
//	 * <BR><I>seed</I> = Random seed
//	 * <BR><I>N</I> = Number of data points
//	 */
//	public static void main
//		(String[] args)
//		throws Exception
//		{
//		if (args.length != 2) usage();
//		long seed = Long.parseLong (args[0]);
//		int N = Integer.parseInt (args[1]);
//		Random prng = Random.getInstance (seed);
//		double[] data = new double [N];
//		for (int i = 0; i < N; ++ i)
//			{
//			data[i] = prng.nextDouble();
//			}
//		Arrays.sort (data);
//		double D = ksTest (data);
//		System.out.println ("D = " + D);
//		System.out.println ("p = " + ksPvalue (N, D));
//		}
//
//	private static void usage()
//		{
//		System.err.println ("Usage: java edu.rit.numeric.Statistics <seed> <N>");
//		System.err.println ("<seed> = Random seed");
//		System.err.println ("<N> = Number of data points");
//		System.exit (1);
//		}

	}
