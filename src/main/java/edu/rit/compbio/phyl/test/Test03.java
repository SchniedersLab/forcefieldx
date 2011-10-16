//*****************************************************************************
//
// File:    Test03.java
// Package: edu.rit.compbio.phyl.test
// Unit:    Class edu.rit.compbio.phyl.test.Test03
//
// This Java source file is copyright (C) 2008 by Alan Kaminsky. All rights
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

package edu.rit.compbio.phyl.test;

import edu.rit.compbio.phyl.DnaSequence;
import edu.rit.compbio.phyl.DnaSequenceTree;

/**
 * Class Test03 is a unit test program for class {@linkplain
 * edu.rit.compbio.phyl.DnaSequenceTree}. It generates all possible rooted
 * bifurcated trees with a given number of tip nodes and prints each tree on the
 * standard output in Newick Standard format.
 * <P>
 * Usage: java edu.rit.compbio.phyl.test.Test03 <I>N</I>
 * <BR><I>N</I> = Number of tip nodes
 *
 * @author  Alan Kaminsky
 * @version 14-Jul-2008
 */
public class Test03
	{

// Prevent construction.

	private Test03()
		{
		}

// Global variables.

	static int N;
	static DnaSequenceTree[] treestack;
	static DnaSequence[] seqstack;

// Main program.

	/**
	 * Unit test main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		if (args.length != 1) usage();
		N = Integer.parseInt (args[0]);
		treestack = new DnaSequenceTree [N];
		for (int i = 0; i < N; ++ i)
			{
			treestack[i] = new DnaSequenceTree (2*i+1);
			}
		seqstack = new DnaSequence [N];
		for (int i = 0; i < N; ++ i)
			{
			seqstack[i] = new DnaSequence (0, 0, nameForLevel (i));
			}
		treestack[0].add (0, seqstack[0]);
		generateTrees (0);
		}

// Hidden operations.

	/**
	 * Generate all trees with one more tip node than the tree at the given
	 * stack level.
	 *
	 * @param  level  Stack level.
	 */
	private static void generateTrees
		(int level)
		{
		DnaSequenceTree tree = treestack[level];
		if (level == N-1)
			{
			System.out.println (tree);
			}
		else
			{
			DnaSequenceTree nexttree = treestack[level+1];
			int len = tree.length();
			for (int i = 0; i < len; ++ i)
				{
				nexttree.copy (tree);
				nexttree.add (i, seqstack[level+1]);
				generateTrees (level+1);
				}
			}
		}

	/**
	 * Returns the DNA sequence name determined by the given level. Level 0 =
	 * <TT>"A"</TT>, level 1 = <TT>"B"</TT>, and so on.
	 *
	 * @param  level  Level.
	 *
	 * @return  DNA sequence name.
	 */
	private static String nameForLevel
		(int level)
		{
		StringBuilder buf = new StringBuilder();
		int v = level + 1;
		while (v > 0)
			{
			buf.insert (0, v2c[v%27]);
			v /= 27;
			}
		return buf.toString();
		}

	private static final char[] v2c = new char[]
		{' ', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
		 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q',
		 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.compbio.phyl.test.Test03 <N>");
		System.err.println ("<N> = Number of tip nodes");
		System.exit (1);
		}

	}
