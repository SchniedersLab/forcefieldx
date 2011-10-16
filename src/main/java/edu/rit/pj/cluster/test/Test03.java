//******************************************************************************
//
// File:    Test03.java
// Package: edu.rit.pj.cluster.test
// Unit:    Class edu.rit.pj.cluster.test.Test03
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

package edu.rit.pj.cluster.test;

import edu.rit.pj.Comm;

import edu.rit.pj.io.StreamFile;

import edu.rit.util.ByteSequence;

import java.io.File;

/**
 * Class Test03 is a unit test main program for class {@linkplain
 * edu.rit.pj.io.StreamFile}. Each backend process reads the input file named on
 * the command line, then prints the file's contents on the standard output. So
 * with <I>K</I> backend processes, <I>K</I> copies of each line of the input
 * file should be printed.
 * <P>
 * Usage: java -Dpj.np=<I>K</I> edu.rit.pj.cluster.test.Test03 <I>file</I>
 *
 * @author  Alan Kaminsky
 * @version 06-Nov-2006
 */
public class Test03
	{

// Prevent construction.

	private Test03()
		{
		}

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		if (args.length != 1) usage();
		Comm.init (args);
		Comm world = Comm.world();
		int size = world.size();
		int rank = world.rank();
		StreamFile file = new StreamFile (new File (args[0]));
		ByteSequence bs = new ByteSequence (file.getInputStream());
		bs.write (System.out);
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<K> edu.rit.pj.cluster.test.Test03 <file>");
		System.exit (1);
		}

	}
