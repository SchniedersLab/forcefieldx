//******************************************************************************
//
// File:    XYSeriesComplex.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.XYSeriesComplex
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

import java.io.PrintStream;
import java.io.PrintWriter;

/**
 * Class XYSeriesComplex is the abstract base class for a series of (<I>x,y</I>)
 * pairs of values where the <I>x</I> values are real (type <TT>double</TT>) and
 * the <I>y</I> values are complex (type <TT>double</TT>).
 *
 * @author  Alan Kaminsky
 * @version 12-Oct-2007
 */
public abstract class XYSeriesComplex
	{

// Exported constructors.

	/**
	 * Construct a new complex XY series.
	 */
	public XYSeriesComplex()
		{
		}

// Exported operations.

	/**
	 * Returns the number of values in this complex XY series.
	 */
	public abstract int length();

	/**
	 * Determine if this complex XY series is empty.
	 *
	 * @return  True if this series is empty (length = 0), false otherwise.
	 */
	public boolean isEmpty()
		{
		return length() == 0;
		}

	/**
	 * Returns the given <I>x</I> value in this complex XY series.
	 *
	 * @param  i  Index.
	 *
	 * @return  The <I>x</I> value in this series at index <TT>i</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range
	 *     <TT>0</TT> .. <TT>length()-1</TT>.
	 */
	public abstract double x
		(int i);

	/**
	 * Returns the real part of the given complex <I>y</I> value in this complex
	 * XY series.
	 *
	 * @param  i  Index.
	 *
	 * @return  The real part of the complex <I>y</I> value in this series at
	 *          index <TT>i</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range
	 *     <TT>0</TT> .. <TT>length()-1</TT>.
	 */
	public abstract double real
		(int i);

	/**
	 * Returns the imaginary part of the given complex <I>y</I> value in this
	 * complex XY series.
	 *
	 * @param  i  Index.
	 *
	 * @return  The imaginary part of the complex <I>y</I> value in this series
	 *          at index <TT>i</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range
	 *     <TT>0</TT> .. <TT>length()-1</TT>.
	 */
	public abstract double imag
		(int i);

	/**
	 * Returns the magnitude of the given complex <I>y</I> value in this complex
	 * XY series. The magnitude is greater than or equal to 0.
	 *
	 * @param  i  Index.
	 *
	 * @return  The magnitude of the complex <I>y</I> value in this series at
	 *          index <TT>i</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range
	 *     <TT>0</TT> .. <TT>length()-1</TT>.
	 */
	public double magnitude
		(int i)
		{
		double a = real (i);
		double b = imag (i);
		double absa = Math.abs (a);
		double absb = Math.abs (b);
		if (absa == 0.0)
			{
			return absb;
			}
		else if (absb == 0.0)
			{
			return absa;
			}
		else if (absa >= absb)
			{
			double bovera = b/a;
			return absa * Math.sqrt (1.0 + bovera*bovera);
			}
		else
			{
			double aoverb = a/b;
			return absb * Math.sqrt (1.0 + aoverb*aoverb);
			}
		}

	/**
	 * Returns the squared magnitude of the given complex <I>y</I> value in this
	 * complex XY series.
	 *
	 * @param  i  Index.
	 *
	 * @return  The squared magnitude of the complex <I>y</I> value in this
	 *          series at index <TT>i</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range
	 *     <TT>0</TT> .. <TT>length()-1</TT>.
	 */
	public double squaredMagnitude
		(int i)
		{
		double a = real (i);
		double b = imag (i);
		double absa = Math.abs (a);
		double absb = Math.abs (b);
		if (absa == 0.0)
			{
			return absb*absb;
			}
		else if (absb == 0.0)
			{
			return absa*absa;
			}
		else if (absa >= absb)
			{
			double bovera = b/a;
			return absa * absa * (1.0 + bovera*bovera);
			}
		else
			{
			double aoverb = a/b;
			return absb * absb * (1.0 + aoverb*aoverb);
			}
		}

	/**
	 * Returns the phase of the given complex <I>y</I> value in this complex XY
	 * series. The phase is in the range -pi to +pi.
	 *
	 * @param  i  Index.
	 *
	 * @return  The phase of the complex <I>y</I> value in this series at index
	 *          <TT>i</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range
	 *     <TT>0</TT> .. <TT>length()-1</TT>.
	 */
	public double phase
		(int i)
		{
		return Math.atan2 (imag (i), real (i));
		}

	/**
	 * Returns a new series consisting of the real parts of the complex <I>y</I>
	 * values in this complex XY series.
	 *
	 * @return  Series of real parts.
	 */
	public Series realSeries()
		{
		final XYSeriesComplex outer = this;
		return new Series()
			{
			public int length()
				{
				return outer.length();
				}
			public double x (int i)
				{
				return outer.real (i);
				}
			};
		}

	/**
	 * Returns a new series consisting of the imaginary parts of the complex
	 * <I>y</I> values in this complex XY series.
	 *
	 * @return  Series of imaginary parts.
	 */
	public Series imagSeries()
		{
		final XYSeriesComplex outer = this;
		return new Series()
			{
			public int length()
				{
				return outer.length();
				}
			public double x (int i)
				{
				return outer.imag (i);
				}
			};
		}

	/**
	 * Returns a new series consisting of the magnitudes of the complex <I>y</I>
	 * values in this complex XY series. Each magnitude is greater than or equal
	 * to 0.
	 *
	 * @return  Series of magnitudes.
	 */
	public Series magnitudeSeries()
		{
		final XYSeriesComplex outer = this;
		return new Series()
			{
			public int length()
				{
				return outer.length();
				}
			public double x (int i)
				{
				return outer.magnitude (i);
				}
			};
		}

	/**
	 * Returns a new series consisting of the squared magnitudes of the complex
	 * <I>y</I> values in this complex XY series.
	 *
	 * @return  Series of squared magnitudes.
	 */
	public Series squaredMagnitudeSeries()
		{
		final XYSeriesComplex outer = this;
		return new Series()
			{
			public int length()
				{
				return outer.length();
				}
			public double x (int i)
				{
				return outer.squaredMagnitude (i);
				}
			};
		}

	/**
	 * Returns a new series consisting of the phases of the complex <I>y</I>
	 * values in this complex XY series. Each phase is in the range -pi to +pi.
	 *
	 * @return  Series of phases.
	 */
	public Series phaseSeries()
		{
		final XYSeriesComplex outer = this;
		return new Series()
			{
			public int length()
				{
				return outer.length();
				}
			public double x (int i)
				{
				return outer.phase (i);
				}
			};
		}

	/**
	 * Returns a new XY series consisting of the <I>x</I> values and the real
	 * parts of the complex <I>y</I> values in this complex XY series.
	 *
	 * @return  XY series of real parts.
	 */
	public XYSeries realXYSeries()
		{
		final XYSeriesComplex outer = this;
		return new XYSeries()
			{
			public int length()
				{
				return outer.length();
				}
			public double x (int i)
				{
				return outer.x (i);
				}
			public double y (int i)
				{
				return outer.real (i);
				}
			};
		}

	/**
	 * Returns a new XY series consisting of the <I>x</I> values and the
	 * imaginary parts of the complex <I>y</I> values in this complex XY series.
	 *
	 * @return  XY series of imaginary parts.
	 */
	public XYSeries imagXYSeries()
		{
		final XYSeriesComplex outer = this;
		return new XYSeries()
			{
			public int length()
				{
				return outer.length();
				}
			public double x (int i)
				{
				return outer.x (i);
				}
			public double y (int i)
				{
				return outer.imag (i);
				}
			};
		}

	/**
	 * Returns a new XY series consisting of the <I>x</I> values and the
	 * magnitudes of the complex <I>y</I> values in this complex XY series. Each
	 * magnitude is greater than or equal to 0.
	 *
	 * @return  XY series of magnitudes.
	 */
	public XYSeries magnitudeXYSeries()
		{
		final XYSeriesComplex outer = this;
		return new XYSeries()
			{
			public int length()
				{
				return outer.length();
				}
			public double x (int i)
				{
				return outer.x (i);
				}
			public double y (int i)
				{
				return outer.magnitude (i);
				}
			};
		}

	/**
	 * Returns a new XY series consisting of the <I>x</I> values and the squared
	 * magnitudes of the complex <I>y</I> values in this complex XY series.
	 *
	 * @return  XY series of squared magnitudes.
	 */
	public XYSeries squaredMagnitudeXYSeries()
		{
		final XYSeriesComplex outer = this;
		return new XYSeries()
			{
			public int length()
				{
				return outer.length();
				}
			public double x (int i)
				{
				return outer.x (i);
				}
			public double y (int i)
				{
				return outer.squaredMagnitude (i);
				}
			};
		}

	/**
	 * Returns a new XY series consisting of the <I>x</I> values and the phases
	 * of the complex <I>y</I> values in this complex XY series. Each phase is
	 * in the range -pi to +pi.
	 *
	 * @return  XY series of phases.
	 */
	public XYSeries phaseXYSeries()
		{
		final XYSeriesComplex outer = this;
		return new XYSeries()
			{
			public int length()
				{
				return outer.length();
				}
			public double x (int i)
				{
				return outer.x (i);
				}
			public double y (int i)
				{
				return outer.phase (i);
				}
			};
		}

	/**
	 * Print this complex XY series on the standard output. Each line of output
	 * consists of the index, the <I>x</I> value, the <I>y</I> value's real
	 * part, the <I>y</I> value's imaginary part, the <I>y</I> value's
	 * magnitude, and the <I>y</I> value's phase, separated by tabs.
	 */
	public void print()
		{
		print (System.out);
		}

	/**
	 * Print this complex XY series on the given print stream. Each line of
	 * output consists of the index, the <I>x</I> value, the <I>y</I> value's
	 * real part, the <I>y</I> value's imaginary part, the <I>y</I> value's
	 * magnitude, and the <I>y</I> value's phase, separated by tabs.
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
			theStream.print (real(i));
			theStream.print ('\t');
			theStream.print (imag(i));
			theStream.print ('\t');
			theStream.print (magnitude(i));
			theStream.print ('\t');
			theStream.println (phase(i));
			}
		}

	/**
	 * Print this complex XY series on the given print writer. Each line of
	 * output consists of the index, the <I>x</I> value, the <I>y</I> value's
	 * real part, the <I>y</I> value's imaginary part, the <I>y</I> value's
	 * magnitude, and the <I>y</I> value's phase, separated by tabs.
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
			theWriter.print (real(i));
			theWriter.print ('\t');
			theWriter.print (imag(i));
			theWriter.print ('\t');
			theWriter.print (magnitude(i));
			theWriter.print ('\t');
			theWriter.println (phase(i));
			}
		}

	}
