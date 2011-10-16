//******************************************************************************
//
// File:    Stdin.java
// Package: edu.rit.pj.job
// Unit:    Class edu.rit.pj.job.Stdin
//
// This Java source file is copyright (C) 2010 by Alan Kaminsky. All rights
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

package edu.rit.pj.job;

/**
 * Enum Stdin enumerates how the standard input stream of a {@linkplain Job}
 * is redirected.
 *
 * @author  Alan Kaminsky
 * @version 07-Oct-2010
 */
enum Stdin
	{

	/**
	 * The standard input is not redirected.
	 */
	NONE,

	/**
	 * The standard input is read from a file. It is an error if the file does
	 * not exist.
	 */
	FILE,

	}
