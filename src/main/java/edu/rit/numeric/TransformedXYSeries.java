//******************************************************************************
//
// File:    TransformedXYSeries.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.TransformedXYSeries
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
 * Class TransformedXYSeries provides an {@linkplain XYSeries} that is formed by
 * transforming the values in an underlying {@linkplain XYSeries}. The
 * underlying series consists of a list of pairs
 * (<I>x</I><SUB><I>i</I></SUB>,&nbsp;<I>y</I><SUB><I>i</I></SUB>). The transformed
 * series consists of a list of pairs
 * (<I>f</I><SUB>1</SUB>(<I>x</I><SUB><I>i</I></SUB>),&nbsp;<I>f</I><SUB>2</SUB>(<I>y</I><SUB><I>i</I></SUB>)),
 * where each of <I>f</I><SUB>1</SUB> and <I>f</I><SUB>2</SUB> is computed by a
 * {@link Function </CODE>Function<CODE>}.
 *
 * @author  Alan Kaminsky
 * @version 27-Jul-2007
 */
public class TransformedXYSeries
	extends XYSeries
	{

// Hidden data members.

	private XYSeries mySeries;
	private Function myXFunction;
	private Function myYFunction;

// Exported constructors.

	/**
	 * Construct a new transformed X-Y series.
	 *
	 * @param  theSeries     Underlying X-Y series.
	 * @param  theXFunction  Function for transforming the X values, or null not
	 *                       to transform the X values.
	 * @param  theYFunction  Function for transforming the Y values, or null not
	 *                       to transform the Y values.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSeries</TT> is null.
	 */
	public TransformedXYSeries
		(XYSeries theSeries,
		 Function theXFunction,
		 Function theYFunction)
		{
		if (theSeries == null)
			{
			throw new NullPointerException
				("TransformedXYSeries(): theSeries is null");
			}
		mySeries = theSeries;
		myXFunction = theXFunction;
		myYFunction = theYFunction;
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

	}
