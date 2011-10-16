//******************************************************************************
//
// File:    TransformedXYZSeries.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.TransformedXYZSeries
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
 * Class TransformedXYZSeries provides an {@linkplain XYZSeries} that is formed
 * by transforming the values in an underlying {@linkplain XYZSeries}. The
 * underlying series consists of a list of triples
 * (<I>x</I><SUB><I>i</I></SUB>,&nbsp;<I>y</I><SUB><I>i</I></SUB>,&nbsp;<I>z</I><SUB><I>i</I></SUB>).
 * The transformed series consists of a list of triples
 * (<I>f</I><SUB>1</SUB>(<I>x</I><SUB><I>i</I></SUB>),&nbsp;<I>f</I><SUB>2</SUB>(<I>y</I><SUB><I>i</I></SUB>),&nbsp;<I>f</I><SUB>3</SUB>(<I>z</I><SUB><I>i</I></SUB>)),
 * where each of <I>f</I><SUB>1</SUB>, <I>f</I><SUB>2</SUB>, and
 * <I>f</I><SUB>3</SUB> is computed by a {@link Function </CODE>Function<CODE>}.
 *
 * @author  Alan Kaminsky
 * @version 27-Jul-2007
 */
public class TransformedXYZSeries
	extends XYZSeries
	{

// Hidden data members.

	private XYZSeries mySeries;
	private Function myXFunction;
	private Function myYFunction;
	private Function myZFunction;

// Exported constructors.

	/**
	 * Construct a new transformed X-Y-Z series.
	 *
	 * @param  theSeries     Underlying X-Y-Z series.
	 * @param  theXFunction  Function for transforming the X values, or null not
	 *                       to transform the X values.
	 * @param  theYFunction  Function for transforming the Y values, or null not
	 *                       to transform the Y values.
	 * @param  theZFunction  Function for transforming the Z values, or null not
	 *                       to transform the Z values.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSeries</TT> is null.
	 */
	public TransformedXYZSeries
		(XYZSeries theSeries,
		 Function theXFunction,
		 Function theYFunction,
		 Function theZFunction)
		{
		if (theSeries == null)
			{
			throw new NullPointerException
				("TransformedXYZSeries(): theSeries is null");
			}
		mySeries = theSeries;
		myXFunction = theXFunction;
		myYFunction = theYFunction;
		myZFunction = theZFunction;
		}

// Exported operations.

	/**
	 * Returns the number of values in this series.
	 *
	 * @return  Length.
	 */
	public int length()
		{
		return mySeries.length();
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
		double xval = mySeries.x(i);
		if (myXFunction != null) xval = myXFunction.f (xval);
		return xval;
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
		double yval = mySeries.y(i);
		if (myYFunction != null) yval = myYFunction.f (yval);
		return yval;
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
		double zval = mySeries.z(i);
		if (myZFunction != null) zval = myZFunction.f (zval);
		return zval;
		}

	}
