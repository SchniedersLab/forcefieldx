//******************************************************************************
//
// File:    ListSeries.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.ListSeries
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
 * Class ListSeries provides a series of real values (type <TT>double</TT>);
 * the series is variable-length.
 *
 * @author  Alan Kaminsky
 * @version 22-Jul-2007
 */
public class ListSeries
	extends Series
	implements Externalizable
	{

// Hidden data members.

	private static final long serialVersionUID = 1521857061875435370L;

	private static final int INCR = 1024;

	private double[] xArray = new double [INCR];
	private int myLength = 0;

// Exported constructors.

	/**
	 * Construct a new zero-length series.
	 */
	public ListSeries()
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
				("ListSeries.x(): Index = " + i + " out of bounds");
			}
		return xArray[i];
		}

	/**
	 * Clear this series.
	 *
	 * @return  This series.
	 */
	public ListSeries clear()
		{
		myLength = 0;
		return this;
		}

	/**
	 * Add the given data to this series.
	 *
	 * @param  x  X value.
	 *
	 * @return  This series.
	 */
	public ListSeries add
		(double x)
		{
		allocate (1, INCR);
		xArray[myLength] = x;
		++ myLength;
		return this;
		}

	/**
	 * Add the given data array to this series. Each element of <TT>x</TT> is
	 * added to this series as though by the <TT>add(double)</TT> method.
	 *
	 * @param  x  Array of X values.
	 *
	 * @return  This series.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>x</TT> is null.
	 */
	public ListSeries add
		(double[] x)
		{
		return add (x, 0, x.length);
		}

	/**
	 * Add a portion of the given data array to this series. Each element of
	 * <TT>x</TT> from index <TT>off</TT> through index <TT>off+len-1</TT>
	 * inclusive is added to this series as though by the <TT>add(double)</TT>
	 * method.
	 *
	 * @param  x    Array of X values.
	 * @param  off  Index of first element to add.
	 * @param  len  Number of elements to add.
	 *
	 * @return  This series.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>x</TT> is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>off</TT> &lt; 0, <TT>len</TT>
	 *     &lt; 0, or <TT>off+len</TT> &gt; <TT>x.length</TT>.
	 */
	public ListSeries add
		(double[] x,
		 int off,
		 int len)
		{
		if (off < 0 || len < 0 || off + len > x.length)
			{
			throw new IndexOutOfBoundsException();
			}
		allocate (len, len);
		System.arraycopy (x, off, xArray, myLength, len);
		myLength += len;
		return this;
		}

	/**
	 * Add the given {@linkplain Series} to this series. Each X value in the
	 * given series is added to this series as though by the
	 * <TT>add(double)</TT> method.
	 *
	 * @param  theSeries  Series.
	 *
	 * @return  This series.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSeries</TT> is null.
	 */
	public ListSeries add
		(Series theSeries)
		{
		return add (theSeries, 0, theSeries.length());
		}

	/**
	 * Add a portion of the given {@linkplain Series} to this series. Each X
	 * value in the given series from index <TT>off</TT> through index
	 * <TT>off+len-1</TT> inclusive is added to this series as though by the
	 * <TT>add(double)</TT> method.
	 *
	 * @param  theSeries  Series.
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
	public ListSeries add
		(Series theSeries,
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
			++ myLength;
			++ off;
			-- len;
			}
		return this;
		}

	/**
	 * Add data read from the given scanner to this series. Double values are
	 * read from the scanner until there are no more. Each double value is added
	 * to this series as though by the <TT>add(double)</TT> method. The
	 * <TT>add()</TT> method does <I>not</I> close the scanner when finished.
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
	public ListSeries add
		(Scanner scanner)
		{
		while (scanner.hasNextDouble())
			{
			add (scanner.nextDouble());
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
		myLength = n;
		for (int i = 0; i < n; ++ i)
			{
			xArray[i] = in.readDouble();
			}
		}

// Hidden operations.

	/**
	 * Increase the storage allocation for this list series if necessary.
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
			}
		}

	}
