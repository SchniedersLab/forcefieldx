//******************************************************************************
//
// File:    RangeException.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.RangeException
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
 * Class RangeException is an unchecked runtime exception thrown if a function's
 * value is outside the range of the return type.
 *
 * @author  Alan Kaminsky
 * @version 06-Jul-2007
 */
public class RangeException
	extends NumericRuntimeException
	{

// Exported constructors.

	/**
	 * Construct a new range exception with no detail message and no cause.
	 */
	public RangeException()
		{
		super();
		}

	/**
	 * Construct a new range exception with the given detail message and no
	 * cause.
	 *
	 * @param  message  Detail message.
	 */
	public RangeException
		(String message)
		{
		super (message);
		}

	/**
	 * Construct a new range exception with the given cause and the default
	 * detail message.
	 *
	 * @param  cause  Cause.
	 */
	public RangeException
		(Throwable cause)
		{
		super (cause);
		}

	/**
	 * Construct a new range exception with the given detail message and the
	 * given cause.
	 *
	 * @param  message  Detail message.
	 * @param  cause    Cause.
	 */
	public RangeException
		(String message,
		 Throwable cause)
		{
		super (message, cause);
		}

	}
