//******************************************************************************
//
// File:    AggregateXYZSeries.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.AggregateXYZSeries
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

/**
 * Class AggregateXYZSeries provides an {@linkplain XYZSeries} formed by
 * aggregating three {@linkplain Series} together.
 *
 * @author  Alan Kaminsky
 * @version 11-Jun-2007
 */
public class AggregateXYZSeries
	extends XYZSeries
	{

// Hidden data members.

	private Series myXSeries;
	private Series myYSeries;
	private Series myZSeries;

// Exported constructors.

	/**
	 * Construct a new aggregate XYZ series. It is assumed that the <TT>x</TT>,
	 * <TT>y</TT>, and <TT>z</TT> series are the same length.
	 * <P>
	 * <I>Note:</I> This series object stores <I>references</I> to <TT>x</TT>,
	 * <TT>y</TT>, and <TT>z</TT>. Changing the contents of <TT>x</TT>,
	 * <TT>y</TT>, or <TT>z</TT> will change the contents of this series.
	 *
	 * @param  x  Series of X values.
	 * @param  y  Series of Y values.
	 * @param  z  Series of Z values.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>x</TT>, <TT>y</TT>, or <TT>z</TT>
	 *     is null.
	 */
	public AggregateXYZSeries
		(Series x,
		 Series y,
		 Series z)
		{
		if (x == null || y == null || z == null)
			{
			throw new NullPointerException();
			}
		myXSeries = x;
		myYSeries = y;
		myZSeries = z;
		}

// Exported operations.

	/**
	 * Returns the number of values in this series.
	 *
	 * @return  Length.
	 */
	public int length()
		{
		return myXSeries.length();
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
		return myXSeries.x (i);
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
		return myYSeries.x (i);
		}

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
	public double z
		(int i)
		{
		return myZSeries.x (i);
		}

	}
