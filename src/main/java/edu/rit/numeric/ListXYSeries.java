//******************************************************************************
//
// File:    ListXYSeries.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.ListXYSeries
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

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

import java.util.Scanner;

/**
 * Class ListXYSeries provides a series of (<I>x,y</I>) pairs of real values
 * (type <TT>double</TT>); the series is variable-length.
 *
 * @author  Alan Kaminsky
 * @version 22-Jul-2007
 */
public class ListXYSeries
	extends XYSeries
	implements Externalizable
	{

// Hidden data members.

	private static final long serialVersionUID = 3623908540262724935L;

	private static final int INCR = 16;

	private double[] xArray = new double [INCR];
	private double[] yArray = new double [INCR];
	private int myLength = 0;

// Exported constructors.

	/**
	 * Construct a new zero-length XY series.
	 */
	public ListXYSeries()
		{
		}

// Exported operations.

	/**
	 * Returns the number of values in this series.
	 *
	 * @return  Length.
	 */
	public int length()
		{
		return myLength;
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
	public double x
		(int i)
		{
		if (0 > i || i >= myLength)
			{
			throw new ArrayIndexOutOfBoundsException
				("ListXYSeries.x(): Index = " + i + " out of bounds");
			}
		return xArray[i];
		}

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
	public double y
		(int i)
		{
		if (0 > i || i >= myLength)
			{
			throw new ArrayIndexOutOfBoundsException
				("ListXYSeries.y(): Index = " + i + " out of bounds");
			}
		return yArray[i];
		}

	/**
	 * Clear this series.
	 *
	 * @return  This series.
	 */
	public ListXYSeries clear()
		{
		myLength = 0;
		return this;
		}

	/**
	 * Add the given data to this series.
	 *
	 * @param  x  X value.
	 * @param  y  Y value.
	 *
	 * @return  This series.
	 */
	public ListXYSeries add
		(double x,
		 double y)
		{
		allocate (1, INCR);
		xArray[myLength] = x;
		yArray[myLength] = y;
		++ myLength;
		return this;
		}

	/**
	 * Add the given data to this series. The data values are grouped into
	 * pairs, where the first data value is the X value and the second data
	 * value is the Y value. Each pair of data values is added to this series as
	 * though by the <TT>add(double,double)</TT> method. It is assumed that
	 * there are an even number of data values.
	 *
	 * @param  data  Data values.
	 *
	 * @return  This series.
	 */
	public ListXYSeries add
		(double... data)
		{
		int len = data.length / 2;
		allocate (len, len);
		for (int i = 0; i < len; ++ i)
			{
			xArray[myLength] = data[2*i];
			yArray[myLength] = data[2*i+1];
			++ myLength;
			}
		return this;
		}

	/**
	 * Add the given data arrays to this series. Each element of <TT>x</TT> and
	 * <TT>y</TT> is added to this series as though by the
	 * <TT>add(double,double)</TT> method. It is assumed that <TT>x</TT> and
	 * <TT>y</TT> are the same length.
	 *
	 * @param  x  Array of X values.
	 * @param  y  Array of Y values.
	 *
	 * @return  This series.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>x</TT> or <TT>y</TT> is null.
	 */
	public ListXYSeries add
		(double[] x,
		 double[] y)
		{
		return add (x, y, 0, x.length);
		}

	/**
	 * Add a portion of the given data arrays to this series. Each element of
	 * <TT>x</TT> and <TT>y</TT> from index <TT>off</TT> through index
	 * <TT>off+len-1</TT> inclusive is added to this series as though by the
	 * <TT>add(double,double)</TT> method.
	 *
	 * @param  x    Array of X values.
	 * @param  y    Array of Y values.
	 * @param  off  Index of first element to add.
	 * @param  len  Number of elements to add.
	 *
	 * @return  This series.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>x</TT> or <TT>y</TT> is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>off</TT> &lt; 0, <TT>len</TT>
	 *     &lt; 0, <TT>off+len</TT> &gt; <TT>x.length</TT>, or <TT>off+len</TT>
	 *     &gt; <TT>y.length</TT>.
	 */
	public ListXYSeries add
		(double[] x,
		 double[] y,
		 int off,
		 int len)
		{
		if (off < 0 || len < 0 || off + len > x.length || off + len > y.length)
			{
			throw new IndexOutOfBoundsException();
			}
		allocate (len, len);
		System.arraycopy (x, off, xArray, myLength, len);
		System.arraycopy (y, off, yArray, myLength, len);
		myLength += len;
		return this;
		}

	/**
	 * Add the given {@linkplain XYSeries} to this series. Each corresponding X
	 * and Y value in the given X-Y series is added to this series as though by
	 * the <TT>add(double,double)</TT> method.
	 *
	 * @param  theSeries  X-Y series.
	 *
	 * @return  This series.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSeries</TT> is null.
	 */
	public ListXYSeries add
		(XYSeries theSeries)
		{
		return add (theSeries, 0, theSeries.length());
		}

	/**
	 * Add a portion of the given {@linkplain XYSeries} to this series. Each
	 * corresponding X and Y value in the given X-Y series from index
	 * <TT>off</TT> through index <TT>off+len-1</TT> inclusive is added to this
	 * series as though by the <TT>add(double,double)</TT> method.
	 *
	 * @param  theSeries  X-Y series.
	 * @param  off        Index of first element to add.
	 * @param  len        Number of elements to add.
	 *
	 * @return  This series.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSeries</TT> is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>off</TT> &lt; 0, <TT>len</TT>
	 *     &lt; 0, or <TT>off+len</TT> &gt; <TT>theSeries.length()</TT>.
	 */
	public ListXYSeries add
		(XYSeries theSeries,
		 int off,
		 int len)
		{
		if (off < 0 || len < 0 || off + len > theSeries.length())
			{
			throw new IndexOutOfBoundsException();
			}
		allocate (len, len);
		while (len > 0)
			{
			xArray[myLength] = theSeries.x (off);
			yArray[myLength] = theSeries.y (off);
			++ myLength;
			++ off;
			-- len;
			}
		return this;
		}

	/**
	 * Add data read from the given scanner to this series. Double values are
	 * read from the scanner until there are no more. Each pair of consecutive
	 * double values is added to this series as though by the
	 * <TT>add(double,double)</TT> method. Any leftover double value is
	 * discarded. The <TT>add()</TT> method does <I>not</I> close the scanner
	 * when finished.
	 * <P>
	 * To read data from a file, pass a scanner constructed on top of the file;
	 * to read data from an input stream, pass a scanner constructed on top of
	 * the input stream; and so on.
	 *
	 * @param  scanner  Scanner.
	 *
	 * @return  This series.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>scanner</TT> is null.
	 */
	public ListXYSeries add
		(Scanner scanner)
		{
		while (scanner.hasNextDouble())
			{
			double x = scanner.nextDouble();
			if (scanner.hasNextDouble())
				{
				double y = scanner.nextDouble();
				add (x, y);
				}
			}
		return this;
		}

	/**
	 * Write this series to the given object output stream.
	 *
	 * @param  out  Object output stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void writeExternal
		(ObjectOutput out)
		throws IOException
		{
		int n = myLength;
		out.writeInt (n);
		for (int i = 0; i < n; ++ i)
			{
			out.writeDouble (xArray[i]);
			out.writeDouble (yArray[i]);
			}
		}

	/**
	 * Read this series from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException
		{
		int n = in.readInt();
		xArray = new double [n];
		yArray = new double [n];
		myLength = n;
		for (int i = 0; i < n; ++ i)
			{
			xArray[i] = in.readDouble();
			yArray[i] = in.readDouble();
			}
		}

// Hidden operations.

	/**
	 * Increase the storage allocation for this list XY series if necessary.
	 *
	 * @param  len   Number of items to be added.
	 * @param  incr  Amount to increase allocation.
	 */
	private void allocate
		(int len,
		 int incr)
		{
		if (myLength + len > xArray.length)
			{
			double[] newxArray = new double [myLength + incr];
			System.arraycopy (xArray, 0, newxArray, 0, myLength);
			xArray = newxArray;
			double[] newyArray = new double [myLength + incr];
			System.arraycopy (yArray, 0, newyArray, 0, myLength);
			yArray = newyArray;
			}
		}

	}
