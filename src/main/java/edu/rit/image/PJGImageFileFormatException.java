//******************************************************************************
//
// File:    PJGImageFileFormatException.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.PJGImageFileFormatException
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

package edu.rit.image;

import java.io.IOException;

/**
 * Class PJGImageFileFormatException is an IOException thrown if a PJG image
 * file's contents are not formatted correctly. For further information about
 * the PJG format, see class {@linkplain PJGImage}.
 *
 * @author  Alan Kaminsky
 * @version 01-Nov-2007
 */
public class PJGImageFileFormatException
	extends IOException
	{

// Exported constructors.

	/**
	 * Construct a new PJG image file format exception with no detail message
	 * and no cause.
	 */
	public PJGImageFileFormatException()
		{
		super();
		}

	/**
	 * Construct a new PJG image file format exception with the given detail
	 * message and no cause.
	 *
	 * @param  message  Detail message.
	 */
	public PJGImageFileFormatException
		(String message)
		{
		super (message);
		}

	/**
	 * Construct a new PJG image file format exception with the given cause and
	 * the default detail message.
	 *
	 * @param  cause  Cause.
	 */
	public PJGImageFileFormatException
		(Throwable cause)
		{
		super();
		initCause (cause);
		}

	/**
	 * Construct a new PJG image file format exception with the given detail
	 * message and the given cause.
	 *
	 * @param  message  Detail message.
	 * @param  cause    Cause.
	 */
	public PJGImageFileFormatException
		(String message,
		 Throwable cause)
		{
		super (message);
		initCause (cause);
		}

	}
