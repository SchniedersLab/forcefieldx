//******************************************************************************
//
// File:    RobustFit.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.RobustFit
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

import edu.rit.util.Random;

import static java.lang.Math.*;

import java.util.Arrays;
import java.util.Comparator;

/**
 * Class RobustFit uses a robust estimation procedure to fit a series of
 * (<I>x,y</I>) data points to a model. The data series is an instance of class
 * {@linkplain XYSeries}. The model is represented by a {@linkplain
 * ParameterizedFunction} that computes the <I>y</I> value, given an <I>x</I>
 * value. The model also has <I>parameters</I>.
 * <P>
 * Given a data series, a model function, and an initial guess for the parameter
 * values, class RobustFit's <TT>fit()</TT> method finds parameter values that
 * minimize the following <I>metric:</I>
 * <P>
 * &nbsp;&nbsp;&nbsp;&nbsp;<B>&Sigma;</B><SUB><I>i</I></SUB>
 * <I>&rho;</I> (<I>y</I><SUB><I>i</I></SUB> &minus;
 * <I>f</I> (<I>x</I><SUB><I>i</I></SUB>, parameters))
 * <P>
 * where <I>f</I> is the model function and <I>&rho;</I> is one of these metric
 * functions:
 * <UL>
 * <P><LI>
 * Normal: <I>&rho;</I> (<I>z</I>) = <I>z</I><SUP>2</SUP>/2
 * <P><LI>
 * Exponential: <I>&rho;</I> (<I>z</I>) = |<I>z</I>|
 * <P><LI>
 * Cauchy (default): <I>&rho;</I> (<I>z</I>) = log (1 + <I>z</I><SUP>2</SUP>/2)
 * </UL>
 * <P>
 * In other words, the <TT>fit()</TT> method <I>fits</I> the model to the data
 * by adjusting the parameters to minimize the metric.
 * <P>
 * The metric function is the negative logarithm of the probability distribution
 * of the errors in the <I>y</I> values. The above metric functions correspond
 * to normal, two-sided exponential, and Cauchy error distributions.
 * <P>
 * The metric functions differ in how they treat <I>outliers,</I> i.e., data
 * points that deviate from the model. The normal metric function gives
 * increasing weights to points with increasing deviations. However, because of
 * the increasing weights, outlier points may skew the fit (hence, this is not
 * really a "robust" metric function). The exponential metric function gives
 * equal weights to all points, regardless of deviation. This reduces the
 * influence of outliers on the fit, yielding a more robust fit. With the
 * Cauchy metric function, the weights first increase, then decrease as the
 * deviations increase. This reduces the influence of outliers even further.
 * <P>
 * The <TT>fit()</TT> method uses class {@linkplain
 * MDMinimizationDownhillSimplex} to find the parameter values that minimize the
 * metric. The inputs to and outputs from the <TT>fit()</TT> method are stored
 * in fields of an instance of class RobustFit.
 * <P>
 * The <TT>fitWithDistribution()</TT> method uses the <I>bootstrapping</I>
 * technique to determine the distribution of the model parameters, which
 * depends on the error distribution of the data points. Bootstrapping performs
 * multiple iterations of the model fitting procedure. On each iteration, a
 * trial data set the same size as the original data set is created by sampling
 * the original data points with replacement, and model parameters for the trial
 * data set are computed. The <TT>fitWithDistribution()</TT> method outputs a
 * series of the parameter values found at each iteration; the confidence region
 * for the parameters; and the goodness-of-fit <I>p</I>-value.
 *
 * @author  Alan Kaminsky
 * @version 22-Oct-2010
 */
public class RobustFit
	{

// Exported data members.

	/**
	 * The model function. When <TT>model.f()</TT> is called, the <TT>x</TT>
	 * argument is <I>x</I><SUB><I>i</I></SUB>, the <I>x</I> value of a data
	 * point; the <TT>p</TT> argument contains the model parameters; and the
	 * return value is
	 * <I>f</I>&nbsp;(<I>x</I><SUB><I>i</I></SUB>,&nbsp;parameters).
	 */
	public final ParameterizedFunction model;

	/**
	 * The number of parameters in the model, <I>M</I>.
	 */
	public final int M;

	/**
	 * The metric function. By default, this is <TT>CAUCHY</TT>. It can instead
	 * be set to <TT>NORMAL</TT>, <TT>EXPONENTIAL</TT>, or some other metric
	 * function.
	 */
	public Function metric = CAUCHY;

	/**
	 * The model parameters. On input to the <TT>fit()</TT> and
	 * <TT>fitWithDistribution()</TT> methods, <TT>param</TT> contains the
	 * initial guess for the model parameters. On output from the <TT>fit()</TT>
	 * and <TT>fitWithDistribution()</TT> methods, <TT>param</TT> contains the
	 * fitted parameter values.
	 */
	public final double[] param;

	/**
	 * The data series. It contains the (<I>x,y</I>) data points to be fitted to
	 * the model. It is specified as an argument of the <TT>fit()</TT> and
	 * <TT>fitWithDistribution()</TT> methods.
	 */
	public XYSeries data;

	/**
	 * The metric value. An output of the <TT>fit()</TT> and
	 * <TT>fitWithDistribution()</TT> methods. It is set to the value of the
	 * metric for the model with the fitted parameters stored in <TT>param</TT>.
	 */
	public double metricValue;

	/**
	 * The model parameter distribution. An output of the
	 * <TT>fitWithDistribution()</TT> method. <TT>paramSeries</TT> is a
	 * <I>T</I>-element array, where <I>T</I> is the number of trials. Each
	 * element of <TT>paramSeries</TT> is an <I>M</I>-element array giving the
	 * fitted parameter values for the corresponding trial.
	 */
	public double[][] paramSeries;

	/**
	 * The metric values for the model parameter distribution. An output of the
	 * <TT>fitWithDistribution()</TT> method. <TT>metricSeries</TT> is a
	 * <I>T</I>-element array, where <I>T</I> is the number of trials. Each
	 * element of <TT>metricSeries</TT> gives the value of the metric for the
	 * model with the parameters stored in the corresponding element of
	 * <TT>paramSeries</TT>.
	 */
	public double[] metricSeries;

	/**
	 * The lower bound of the confidence region for the model parameters. An
	 * output of the <TT>fitWithDistribution()</TT> method. The confidence level
	 * is specified as an argument of the <TT>fitWithDistribution()</TT> method;
	 * for example, 0.90 specifies a 90% confidence level. The confidence region
	 * is an <I>M</I>-dimensional rectangular hyperprism centered on the fitted
	 * parameters stored in <TT>param</TT>, such that the given fraction of the
	 * model parameter distribution stored in <TT>paramSeries</TT> falls within
	 * the hyperprism. <TT>confidenceRegionLowerBound</TT> gives the lower bound
	 * of each dimension of the confidence region hyperprism.
	 */
	public double[] confidenceRegionLowerBound;

	/**
	 * The upper bound of the confidence region for the model parameters. An
	 * output of the <TT>fitWithDistribution()</TT> method.
	 * <TT>confidenceRegionUpperBound</TT> gives the upper bound of each
	 * dimension of the confidence region hyperprism.
	 */
	public double[] confidenceRegionUpperBound;

	/**
	 * The goodness-of-fit <I>p</I>-value. An output of the
	 * <TT>fitWithDistribution()</TT> method. This gives the probability that a
	 * metric value greater than or equal to <TT>metricValue</TT> would occur by
	 * chance, even if the model with parameters <TT>params</TT> is correct.
	 */
	public double pValue;

// Exported constants.

	/**
	 * The normal metric function.
	 */
	public static final Function NORMAL = new Function()
		{
		public double f (double z)
			{
			return 0.5*z*z;
			}
		};

	/**
	 * The exponential metric function.
	 */
	public static final Function EXPONENTIAL = new Function()
		{
		public double f (double z)
			{
			return abs(z);
			}
		};

	/**
	 * The Cauchy metric function.
	 */
	public static final Function CAUCHY = new Function()
		{
		public double f (double z)
			{
			return log(1.0 + 0.5*z*z);
			}
		};

// Hidden data members.

	// For minimizing the metric.
	private MDMinimizationDownhillSimplex minimizer;

// Hidden helper classes.

	/**
	 * Class Metric computes the metric to be minimized.
	 *
	 * @author  Alan Kaminsky
	 * @version 06-Oct-2010
	 */
	private class Metric
		implements MDFunction
		{
		public int argumentLength()
			{
			return M; // Number of model parameters
			}

		public double f
			(double[] p) // Model parameters
			{
			double sum = 0.0;
			int N = data.length();
			for (int i = 0; i < N; ++ i)
				{
				sum += metric.f (data.y(i) - model.f (data.x(i), p));
				}
			return sum;
			}
		}

// Exported constructors.

	/**
	 * Construct a new robust fitting object for the given model. The
	 * <TT>model</TT> field is set to the corresponding argument. The <TT>M</TT>
	 * field is set by calling the model function's <TT>parameterLength()</TT>
	 * method. The <TT>param</TT> field is allocated with <I>M</I> elements;
	 * initially, the elements are 0.
	 *
	 * @param  model  Model function.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>model</TT> is null.
	 */
	public RobustFit
		(ParameterizedFunction model)
		{
		if (model == null)
			{
			throw new NullPointerException
				("RobustFit(): model is null");
			}

		this.model = model;
		this.M = model.parameterLength();
		this.param = new double [M];

		minimizer = new MDMinimizationDownhillSimplex (new Metric());
		}

// Exported operations.

	/**
	 * Fit the given data series to the model. The data series is stored in the
	 * <TT>data</TT> field. The model function was specified to the constructor,
	 * and is also stored in the <TT>model</TT> field. On input to the
	 * <TT>fit()</TT> method, <TT>param</TT> contains the initial guess for the
	 * model parameters. On output from the <TT>fit()</TT> method,
	 * <TT>param</TT> contains the fitted parameter values and
	 * <TT>metricValue</TT> contains the value of the metric for the fitted
	 * parameters.
	 * <P>
	 * The <TT>fit()</TT> method uses the downhill simplex technique to find the
	 * model parameters that minimize the metric. This involves initializing the
	 * <I>simplex</I> in an {@linkplain MDMinimizationDownhillSimplex} object.
	 * The <TT>initializeSimplex()</TT> method is called to initialize the
	 * simplex.
	 *
	 * @param  data  Data series.
	 *
	 * @exception  TooManyIterationsException
	 *     (unchecked exception) Thrown if too many iterations occurred without
	 *     finding parameters that minimize the metric function.
	 */
	public void fit
		(XYSeries data)
		{
		this.data = data;
		initializeSimplex (minimizer);
		minimizer.minimize();
		System.arraycopy (minimizer.x[0], 0, param, 0, M);
		metricValue = minimizer.f[0];
		}

	/**
	 * Fit the given data series to the model and compute the distribution of
	 * the model parameters. The data series is stored in the <TT>data</TT>
	 * field. The bootstrapping technique with <I>T</I> trials using the given
	 * pseudorandom number generator is used to compute the distribution. The
	 * given confidence level is used to compute the confidence region; for
	 * example, 0.90 specifies a 90% confidence level. The model function was
	 * specified to the constructor, and is also stored in the <TT>model</TT>
	 * field. On input to the <TT>fitWithDistribution()</TT> method,
	 * <TT>param</TT> contains the initial guess for the model parameters. On
	 * output from the <TT>fit()</TT> method, <TT>param</TT> contains the fitted
	 * parameter values, <TT>metricValue</TT> contains the value of the metric
	 * for the fitted parameters, <TT>paramSeries</TT> contains the series of
	 * fitted parameter values from all the trials, <TT>metricSeries</TT>
	 * contains the metric values from all the trials,
	 * <TT>confidenceRegionLowerBound</TT> and
	 * <TT>confidenceRegionUpperBound</TT> contain the lower and upper bounds of
	 * the confidence region hyperprism, and <TT>pValue</TT> contains the
	 * goodness-of-fit.
	 * <P>
	 * The <TT>fitWithDistribution()</TT> method uses the downhill simplex
	 * technique to find the model parameters that minimize the metric. This
	 * involves initializing the <I>simplex</I> in an {@linkplain
	 * MDMinimizationDownhillSimplex} object. The <TT>initializeSimplex()</TT>
	 * method is called to initialize the simplex.
	 *
	 * @param  data  Data series.
	 * @param  T     Number of trials.
	 * @param  prng  Pseudorandom number generator.
	 * @param  conf  Confidence level, in the range 0.0 .. 1.0.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>conf</TT> is out of bounds.
	 * @exception  TooManyIterationsException
	 *     (unchecked exception) Thrown if too many iterations occurred without
	 *     finding parameters that minimize the metric function.
	 */
	public void fitWithDistribution
		(XYSeries data,
		 int T,
		 Random prng,
		 double conf)
		{
		if (0.0 > conf || conf > 1.0)
			{
			throw new IllegalArgumentException
				("RobustFit.fitWithDistribution(): conf = "+conf+" illegal");
			}

		int N = data.length();

		// Allocate storage for outputs.
		paramSeries = new double [T] [M];
		metricSeries = new double [T];
		confidenceRegionLowerBound = new double [M];
		confidenceRegionUpperBound = new double [M];

		// Compute fitted parameters for original data set.
		this.data = data;
		initializeSimplex (minimizer);
		minimizer.minimize();
		System.arraycopy (minimizer.x[0], 0, param, 0, M);
		metricValue = minimizer.f[0];

		// Perform bootstrap trials.
		ListXYSeries trialData = new ListXYSeries();
		this.data = trialData;
		for (int trial = 0; trial < T; ++ trial)
			{
			// Create trial data set.
			trialData.clear();
			for (int i = 0; i < N; ++ i)
				{
				int j = prng.nextInt (N);
				trialData.add (data.x(j), data.y(j));
				}

			// Perform robust fit on trial data set.
			initializeSimplex (minimizer);
			minimizer.minimize();

			// Record fitted parameters and metric value.
			System.arraycopy (minimizer.x[0], 0, paramSeries[trial], 0, M);
			metricSeries[trial] = minimizer.f[0];
			}
		this.data = data;

		// Compute confidence region: First, compute distance of each
		// bootstrapped parameter set from the fitted parameter set.
		final double[] distance = new double [T];
		for (int trial = 0; trial < T; ++ trial)
			{
			double[] paramSeries_trial = paramSeries[trial];
			double maxDistance = 0.0;
			for (int i = 0; i < M; ++ i)
				{
				maxDistance =
					max (maxDistance, abs (paramSeries_trial[i] - param[i]));
				}
			distance[trial] = maxDistance;
			}

		// Compute confidence region: Second, create an index array for visiting
		// the distance, paramSeries, and metricSeries arrays in ascending order
		// of distance.
		Integer[] index = new Integer [T];
		for (int i = 0; i < T; ++ i) index[i] = i;
		Arrays.sort (index, new Comparator<Integer>()
			{
			public int compare (Integer a, Integer b)
				{
				return
					distance[a] < distance[b] ? -1 :
					distance[a] > distance[b] ? +1 :
					0;
				}
			});

		// Compute confidence region: Third, visit the specified fraction of the
		// bootstrapped parameter sets, in ascending order of distance from the
		// fitted parameter set, and determine confidence region lower and upper
		// bounds.
		for (int i = 0; i < M; ++ i)
			{
			confidenceRegionLowerBound[i] = Double.POSITIVE_INFINITY;
			confidenceRegionUpperBound[i] = Double.NEGATIVE_INFINITY;
			}
		int conf_T = (int)(conf*T + 0.5);
		for (int trial = 0; trial < conf_T; ++ trial)
			{
			double[] paramSeries_trial = paramSeries[index[trial]];
			for (int i = 0; i < M; ++ i)
				{
				confidenceRegionLowerBound[i] =
					min (confidenceRegionLowerBound[i], paramSeries_trial[i]);
				confidenceRegionUpperBound[i] =
					max (confidenceRegionUpperBound[i], paramSeries_trial[i]);
				}
			}

		// Compute p-value.
		pValue = 0.0;
		for (int trial = 0; trial < T; ++ trial)
			{
			if (metricSeries[trial] >= metricValue) pValue += 1.0;
			}
		pValue /= T;
		}

// Hidden operations.

	/**
	 * Initialize the simplex in the given downhill simplex minimizer object.
	 * The simplex points must be set based on the initial guess for the
	 * parameter values stored in the <TT>param</TT> field. For further
	 * information about initializing the simplex, see class {@linkplain
	 * MDMinimizationDownhillSimplex}.
	 * <P>
	 * The default implementation of this method sets the first simplex point to
	 * <TT>param</TT>, sets the second simplex point to <TT>param</TT> except
	 * element 0 is set to perturb(<TT>param[0]</TT>), sets the third simplex
	 * point to <TT>param</TT> except element 1 is set to
	 * perturb(<TT>param[1]</TT>), and so on. perturb(<I>x</I>) = 1.01<I>x</I>
	 * if <I>x</I> &ne; 0; perturb(0) = 0.01. This method can be overridden to
	 * initialize the simplex differently.
	 */
	protected void initializeSimplex
		(MDMinimizationDownhillSimplex minimizer)
		{
		System.arraycopy (param, 0, minimizer.x[0], 0, M);
		for (int i = 1; i <= M; ++ i)
			{
			double[] x_i = minimizer.x[i];
			System.arraycopy (param, 0, x_i, 0, M);
			x_i[i-1] = perturb (x_i[i-1]);
			}
		}

	private static double perturb
		(double x)
		{
		return x == 0.0 ? 0.01 : 1.01*x;
		}

	}
