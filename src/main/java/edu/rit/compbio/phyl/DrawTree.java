//******************************************************************************
//
// File:    DrawTree.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.DrawTree
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

package edu.rit.compbio.phyl;

import edu.rit.draw.Drawing;

import java.io.File;

/**
 * Class DrawTree is a program that draws a phylogenetic tree. The tree is
 * specified on the command line in a subset of Newick Standard format. The
 * drawing is stored in a file specified on the command line as a serialized
 * {@linkplain edu.rit.draw.Drawing} object. The {@linkplain View} program can
 * be used to view the drawing file and save it in an image file in several
 * formats.
 * <P>
 * Usage: java edu.rit.compbio.phyl.DrawTree "<I>tree</I>" <I>file</I> [
 * <I>format</I> ]
 * <BR><I>tree</I> = Tree in Newick Standard format
 * <BR><I>file</I> = Drawing file name
 * <BR><I>format</I> = DecimalFormat string for branch lengths (default:
 * <TT>"0.00"</TT>)
 * <P>
 * <B>Newick Standard format.</B> The subset of Newick Standard format class
 * DrawTree supports is described by this BNF syntax. Nonterminal symbols are
 * in <I>italic</I> font, terminal symbols are in <TT>typewriter</TT> font, |
 * designates alternatives.
 * <P>
 * <I>tree</I> ::= <I>node</I> <TT>;</TT>
 * <BR><I>node</I> ::= <TT>(</TT> <I>childList</I> <TT>)</TT> <I>nameLength</I>
 * <BR><I>childList</I> ::= <I>child</I> | <I>child</I> <TT>,</TT> <I>childList</I>
 * <BR><I>child</I> ::= <I>tip</I> | <I>node</I>
 * <BR><I>tip</I> ::= <I>nameLength</I>
 * <BR><I>nameLength</I> ::= <I>empty</I> | <I>name</I> | <TT>:</TT> <I>length</I> | <I>name</I> <TT>:</TT> <I>length</I>
 * <BR><I>name</I> ::= Any sequence of non-whitespace characters except <TT>(</TT> <TT>)</TT> <TT>,</TT> <TT>:</TT> <TT>;</TT>
 * <BR><I>length</I> ::= Floating point number, as in the <TT>Double.valueOf()</TT> method
 * <P>
 * This subset only supports non-quoted names with no whitespace; however, any
 * underscore character <TT>_</TT> in a name is replaced with a space character.
 * The full Newick Standard format supports quoted names including whitespace.
 * For further information about Newick Standard format, see:
 * <UL>
 * <LI>
 * <A HREF="http://evolution.gs.washington.edu/phylip/newicktree.html" TARGET="_top">http://evolution.gs.washington.edu/phylip/newicktree.html</A>
 * </UL>
 * <P>
 * Here is an example command to draw a tree and store it in the file
 * <TT>"tree.dwg"</TT>:
 * <P>
 * <TT>java edu.rit.compbio.phyl.DrawTree \</TT>
 * <BR><TT>"(Gibbon:10,(Orangutan:8,(Gorilla:6,(Chimp:4,Human:4)))apes)primates;" tree.dwg</TT>
 * <P>
 * Here is the command to view the tree drawing:
 * <P>
 * <TT>java View tree.dwg</TT>
 * <P>
 * Here is the resulting drawing:
 * <P>
 * <IMG SRC="doc-files/tree.png">
 *
 * @author  Alan Kaminsky
 * @version 15-Jul-2008
 */
public class DrawTree
	{

// Prevent construction.

	private DrawTree()
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
		// Parse command line arguments.
		if (args.length < 2 || args.length > 3) usage();
		String treestring = args[0];
		File file = new File (args[1]);
		TreeDrawing artist = new TreeDrawing();
		if (args.length >= 3) artist.setBranchLengthFormat (args[2]);

		// Draw the tree.
		try
			{
			artist.draw (treestring);
			}
		catch (TreeDrawing.SyntaxException exc)
			{
			exc.printSyntaxError (System.err);
			System.exit (1);
			}

		// Write drawing object into output file.
		Drawing.write (file);
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.compbio.phyl.DrawTree \"<tree>\" <file> [<format>]");
		System.err.println ("<tree> = Tree in Newick Standard format");
		System.err.println ("<file> = Drawing file name");
		System.err.println ("<format> = DecimalFormat string for branch lengths (default: \"0.00\")");
		System.exit (1);
		}

	}
