//******************************************************************************
//
// File:    IntegerOp.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.IntegerOp
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

package edu.rit.pj.reduction;

/**
 * Class IntegerOp is the abstract base class for a binary operation on integer
 * values, used to do reduction in a parallel program.
 *
 * @author  Alan Kaminsky
 * @version 24-Nov-2009
 */
public abstract class IntegerOp
	extends Op
	{

// Hidden constructors.

	/**
	 * Construct a new integer binary operation.
	 */
	protected IntegerOp()
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
	public abstract int op
		(int x,
		 int y);

// Exported constants.

	/**
	 * The integer sum binary operation.
	 */
	public static final IntegerOp SUM =
		new IntegerOp()
			{
			public int op
				(int x,
				 int y)
				{
				return x + y;
				}
			};

	/**
	 * The integer product binary operation.
	 */
	public static final IntegerOp PRODUCT =
		new IntegerOp()
			{
			public int op
				(int x,
				 int y)
				{
				return x * y;
				}
			};

	/**
	 * The integer minimum binary operation.
	 */
	public static final IntegerOp MINIMUM =
		new IntegerOp()
			{
			public int op
				(int x,
				 int y)
				{
				return Math.min (x, y);
				}
			};

	/**
	 * The integer maximum binary operation.
	 */
	public static final IntegerOp MAXIMUM =
		new IntegerOp()
			{
			public int op
				(int x,
				 int y)
				{
				return Math.max (x, y);
				}
			};

	/**
	 * The integer bitwise "and" binary operation.
	 */
	public static final IntegerOp AND =
		new IntegerOp()
			{
			public int op
				(int x,
				 int y)
				{
				return x & y;
				}
			};

	/**
	 * The integer bitwise "or" binary operation.
	 */
	public static final IntegerOp OR =
		new IntegerOp()
			{
			public int op
				(int x,
				 int y)
				{
				return x | y;
				}
			};

	/**
	 * The integer bitwise "exclusive or" binary operation.
	 */
	public static final IntegerOp XOR =
		new IntegerOp()
			{
			public int op
				(int x,
				 int y)
				{
				return x ^ y;
				}
			};

	}
