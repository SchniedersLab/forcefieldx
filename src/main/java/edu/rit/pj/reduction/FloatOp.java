//******************************************************************************
//
// File:    FloatOp.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.FloatOp
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

package edu.rit.pj.reduction;

/**
 * Class FloatOp is the abstract base class for a binary operation on float
 * values, used to do reduction in a parallel program.
 *
 * @author  Alan Kaminsky
 * @version 05-Jun-2007
 */
public abstract class FloatOp
	extends Op
	{

// Hidden constructors.

	/**
	 * Construct a new float binary operation.
	 */
	protected FloatOp()
		{
		super();
		}

// Exported operations.

	/**
	 * Perform this binary operation.
	 *
	 * @param  x  First argument.
	 * @param  y  Second argument.
	 *
	 * @return  (<TT>x</TT> <I>op</I> <TT>y</TT>), where <I>op</I> stands for
	 *          this binary operation.
	 */
	public abstract float op
		(float x,
		 float y);

// Exported constants.

	/**
	 * The float sum binary operation.
	 */
	public static final FloatOp SUM =
		new FloatOp()
			{
			public float op
				(float x,
				 float y)
				{
				return x + y;
				}
			};

	/**
	 * The float product binary operation.
	 */
	public static final FloatOp PRODUCT =
		new FloatOp()
			{
			public float op
				(float x,
				 float y)
				{
				return x * y;
				}
			};

	/**
	 * The float minimum binary operation.
	 */
	public static final FloatOp MINIMUM =
		new FloatOp()
			{
			public float op
				(float x,
				 float y)
				{
				return Math.min (x, y);
				}
			};

	/**
	 * The float maximum binary operation.
	 */
	public static final FloatOp MAXIMUM =
		new FloatOp()
			{
			public float op
				(float x,
				 float y)
				{
				return Math.max (x, y);
				}
			};

	}
