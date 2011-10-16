//******************************************************************************
//
// File:    CryptoUtils.java
// Package: edu.rit.crypto
// Unit:    Class edu.rit.crypto.CryptoUtils
//
// This Java source file is copyright (C) 2006 by Alan Kaminsky. All rights
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

package edu.rit.crypto;

/**
 * Class CryptoUtils contains static cryptographic utility methods.
 *
 * @author  Alan Kaminsky
 * @version 05-Apr-2006
 */
public class CryptoUtils
	{

// Prevent construction.

	private CryptoUtils()
		{
		}

// Exported operations.

	/**
	 * Erase the given byte array. Every byte is set to 0.
	 *
	 * @param  theArray  Byte array.
	 */
	public static void erase
		(byte[] theArray)
		{
		int n = theArray == null ? 0 : theArray.length;
		for (int i = 0; i < n; ++ i)
			{
			theArray[i] = (byte) 0;
			}
		}

	/**
	 * Erase the given character array. Every character is set to 0.
	 *
	 * @param  theArray  Character array.
	 */
	public static void erase
		(char[] theArray)
		{
		int n = theArray == null ? 0 : theArray.length;
		for (int i = 0; i < n; ++ i)
			{
			theArray[i] = (char) 0;
			}
		}

	/**
	 * Erase the given integer array. Every integer is set to 0.
	 *
	 * @param  theArray  Integer array.
	 */
	public static void erase
		(int[] theArray)
		{
		int n = theArray == null ? 0 : theArray.length;
		for (int i = 0; i < n; ++ i)
			{
			theArray[i] = 0;
			}
		}

	/**
	 * Erase the given long array. Every long is set to 0.
	 *
	 * @param  theArray  Long array.
	 */
	public static void erase
		(long[] theArray)
		{
		int n = theArray == null ? 0 : theArray.length;
		for (int i = 0; i < n; ++ i)
			{
			theArray[i] = 0L;
			}
		}

	}
