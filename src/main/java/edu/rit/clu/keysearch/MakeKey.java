//******************************************************************************
//
// File:    MakeKey.java
// Package: edu.rit.clu.keysearch
// Unit:    Class edu.rit.clu.keysearch.MakeKey
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

package edu.rit.clu.keysearch;

import edu.rit.util.Hex;

import java.security.SecureRandom;

/**
 * Class MakeKey creates a random 256-bit key for a block cipher. The program
 * gets 32 bytes from the kernel's entropy source using class
 * java.security.SecureRandom's <TT>getSeed()</TT> method, and prints the bytes
 * as a hexadecimal string on the standard output.
 * <P>
 * Usage: java edu.rit.clu.keysearch.MakeKey
 *
 * @author  Alan Kaminsky
 * @version 03-Jul-2007
 */
public class MakeKey
	{

// Prevent construction.

	private MakeKey()
		{
		}

// Main program.

	/**
	 * MakeKey main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		System.out.println (Hex.toString (SecureRandom.getSeed (32)));
		}

	}
