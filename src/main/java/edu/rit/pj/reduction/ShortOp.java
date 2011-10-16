//******************************************************************************
//
// File:    ShortOp.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.ShortOp
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
 * Class ShortOp is the abstract base class for a binary operation on short
 * values, used to do reduction in a parallel program.
 *
 * @author  Alan Kaminsky
 * @version 24-Nov-2009
 */
public abstract class ShortOp
	extends Op
	{

// Hidden constructors.

	/**
	 * Construct a new short binary operation.
	 */
	protected ShortOp()
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
	public abstract short op
		(short x,
		 short y);

// Exported constants.

	/**
	 * The short sum binary operation.
	 */
	public static final ShortOp SUM =
		new ShortOp()
			{
			public short op
				(short x,
				 short y)
				{
				return (short) (x + y);
				}
			};

	/**
	 * The short product binary operation.
	 */
	public static final ShortOp PRODUCT =
		new ShortOp()
			{
			public short op
				(short x,
				 short y)
				{
				return (short) (x * y);
				}
			};

	/**
	 * The short minimum binary operation.
	 */
	public static final ShortOp MINIMUM =
		new ShortOp()
			{
			public short op
				(short x,
				 short y)
				{
				return (short) Math.min (x, y);
				}
			};

	/**
	 * The short maximum binary operation.
	 */
	public static final ShortOp MAXIMUM =
		new ShortOp()
			{
			public short op
				(short x,
				 short y)
				{
				return (short) Math.max (x, y);
				}
			};

	/**
	 * The short bitwise "and" binary operation.
	 */
	public static final ShortOp AND =
		new ShortOp()
			{
			public short op
				(short x,
				 short y)
				{
				return (short)(x & y);
				}
			};

	/**
	 * The short bitwise "or" binary operation.
	 */
	public static final ShortOp OR =
		new ShortOp()
			{
			public short op
				(short x,
				 short y)
				{
				return (short)(x | y);
				}
			};

	/**
	 * The short bitwise "exclusive or" binary operation.
	 */
	public static final ShortOp XOR =
		new ShortOp()
			{
			public short op
				(short x,
				 short y)
				{
				return (short)(x ^ y);
				}
			};

	}
