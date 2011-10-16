//******************************************************************************
//
// File:    ArraySeries.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.ArraySeries
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

package edu.rit.numeric;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class ArraySeries provides a {@linkplain Series} view of values stored in an
 * array. Changing the contents of the array will change the ArraySeries.
 *
 * @author  Alan Kaminsky
 * @version 17-Jun-2008
 */
public class ArraySeries
	extends Series
	implements Externalizable
	{

// Hidden data members.

	private static final long serialVersionUID = 5850541176666797722L;

	private double[] xarray;

// Exported constructors.

	/**
	 * Construct a new uninitialized array series. This constructor is for use
	 * only by object deserialization.
	 */
	public ArraySeries()
		{
		}

	/**
	 * Construct a new array series.
	 *
	 * @param  xarray  Array of X values.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>xarray</TT> is null.
	 */
	public ArraySeries
		(double[] xarray)
		{
		if (xarray == null)
			{
			throw new NullPointerException
				("ArraySeries(): xarray is null");
			}
		this.xarray = xarray;
		}

// Exported operations.

	/**
	 * Returns the number of values in this series.
	 *
	 * @return  Length.
	 */
	public int length()
		{
		return xarray.length;
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
		return xarray[i];
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
		int n = xarray.length;
		out.writeInt (n);
		for (int i = 0; i < n; ++ i)
			{
			out.writeDouble (xarray[i]);
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
		xarray = new double [n];
		for (int i = 0; i < n; ++ i)
			{
			xarray[i] = in.readDouble();
			}
		}

	}
