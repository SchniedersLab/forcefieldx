//******************************************************************************
//
// File:    Constants.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.Constants
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

package edu.rit.mp;

/**
 * Class Constants provides various constants used in the Message Protocol (MP).
 *
 * @author  Alan Kaminsky
 * @version 26-Nov-2007
 */
class Constants
	{

// Prevent construction.

	private Constants()
		{
		}

// Magic number.

	static final int MAGIC_NUMBER = 30144596;

// Item types.

	static final byte TYPE_BOOLEAN   = (byte) 1;
	static final byte TYPE_BYTE      = (byte) 2;
	static final byte TYPE_SHORT     = (byte) 3;
	static final byte TYPE_INTEGER   = (byte) 4;
	static final byte TYPE_LONG      = (byte) 5;
	static final byte TYPE_CHARACTER = (byte) 6;
	static final byte TYPE_FLOAT     = (byte) 7;
	static final byte TYPE_DOUBLE    = (byte) 8;
	static final byte TYPE_OBJECT    = (byte) 9;
	static final byte TYPE_SIGNED_8_BIT_INTEGER    = (byte) 10;
	static final byte TYPE_SIGNED_16_BIT_INTEGER   = (byte) 11;
	static final byte TYPE_UNSIGNED_8_BIT_INTEGER  = (byte) 12;
	static final byte TYPE_UNSIGNED_16_BIT_INTEGER = (byte) 13;

	static final byte TYPE_FIRST = TYPE_BOOLEAN;
	static final byte TYPE_LAST  = TYPE_UNSIGNED_16_BIT_INTEGER;

// Default buffer size in bytes. This is the Ethernet MTU (1500) minus the IP
// header (20) minus the TCP header (20), times 20.

	static final int BUFFER_SIZE = 29200;

	}
