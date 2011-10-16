//******************************************************************************
//
// File:    TreeDrawing.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.TreeDrawing
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

import edu.rit.draw.item.Group;
import edu.rit.draw.item.Line;
import edu.rit.draw.item.Point;
import edu.rit.draw.item.Text;

import java.io.PrintStream;

import java.text.DecimalFormat;

import java.util.LinkedList;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class TreeDrawing provides an object that draws a picture of a phylogenetic
 * tree. The tree is specified either as a {@linkplain DnaSequenceTree} object
 * or as a string in a subset of Newick Standard format.
 * <P>
 * <I>Note:</I> Class TreeDrawing is not multiple thread safe.
 * <P>
 * <B>Newick Standard format.</B> The subset of Newick Standard format class
 * TreeDrawing supports is described by this BNF syntax. Nonterminal symbols are
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
 * Here is an example tree:
 * <P>
 * <TT>"(Gibbon:10,(Orangutan:8,(Gorilla:6,(Chimp:4,Human:4)))apes)primates;"</TT>
 * <P>
 * Here is the resulting drawing:
 * <P>
 * <IMG SRC="doc-files/tree.png">
 *
 * @author  Alan Kaminsky
 * @version 25-Jul-2008
 */
public class TreeDrawing
	{

// Hidden constants.

	// Vertical distance between tip strings.
	private static final double V = 27.0;

	// Length of lines from a node to its descendents.
	private static final double H = 72.0;

	// Gap from line to tip text.
	private static final double GAP = 6.0;

	// Gap from line to branch length text.
	private static final double VGAP = 3.0;

	// Regular expression for parsing tokens in Newick Standard format.
	private static final Pattern pattern =
		Pattern.compile ("\\(|\\)|,|:|;|[^(),:; \\t\\n\\x0B\\f\\r]+");

// Hidden data members.

	// Tree string in Newick Standard format.
	private String treestring;

	// For formatting branch lengths.
	private DecimalFormat format = new DecimalFormat ("0.00");

	// Matcher for parsing tokens.
	private Matcher matcher;

	// Current token and its starting index in <treestring>.
	private String token;
	private int index;

	// Number of tips.
	private int tipCount;

	// Group in which to store the tree drawing.
	private Group group;

// Exported constructors.

	/**
	 * Construct a new tree drawing object.
	 */
	public TreeDrawing()
		{
		}

// Main program.

	/**
	 * Specify the format with which to draw branch lengths. The format string
	 * is used to construct a java.text.DecimalFormat object. If not specified,
	 * the default format string is <TT>"0.00"</TT>.
	 *
	 * @param  format  Format string.
	 */
	public void setBranchLengthFormat
		(String format)
		{
		this.format = new DecimalFormat (format);
		}

	/**
	 * Draw a picture of the given tree. The picture is added to the default
	 * {@linkplain edu.rit.draw.Drawing} object.
	 *
	 * @param  tree  DNA sequence tree.
	 */
	public void draw
		(DnaSequenceTree tree)
		{
		draw (tree, Drawing.defaultDrawing());
		}

	/**
	 * Draw a picture of the given tree. The picture is added to the given
	 * {@linkplain edu.rit.draw.Drawing} object.
	 *
	 * @param  tree  DNA sequence tree.
	 * @param  dwg   Drawing object.
	 */
	public void draw
		(DnaSequenceTree tree,
		 Drawing dwg)
		{
		draw (tree, new Group());
		dwg.add (group);
		}

	/**
	 * Draw a picture of the given tree. The picture is added to the given
	 * {@linkplain edu.rit.draw.item.Group} object.
	 *
	 * @param  tree   DNA sequence tree.
	 * @param  group  Drawing group object.
	 */
	public void draw
		(DnaSequenceTree tree,
		 Group group)
		{
		// Initialize.
		this.tipCount = 0;
		this.group = group;

		// Traverse the tree, building up the drawing as it goes.
		Point root = drawNode (tree, tree.root());

		// Add final line to root.
		group.append (new Line().to(root).hby(-H));
		}

	/**
	 * Draw a picture of the given tree string. The picture is added to the
	 * default {@linkplain edu.rit.draw.Drawing} object.
	 *
	 * @param  tree  Tree as a string in Newick Standard format.
	 *
	 * @exception  SyntaxException
	 *     Thrown if there was a syntax error in <TT>tree</TT>.
	 */
	public void draw
		(String tree)
		throws SyntaxException
		{
		draw (tree, Drawing.defaultDrawing());
		}

	/**
	 * Draw a picture of the given tree string. The picture is added to the
	 * given {@linkplain edu.rit.draw.Drawing} object.
	 *
	 * @param  tree  Tree as a string in Newick Standard format.
	 * @param  dwg   Drawing object.
	 *
	 * @exception  SyntaxException
	 *     Thrown if there was a syntax error in <TT>tree</TT>.
	 */
	public void draw
		(String tree,
		 Drawing dwg)
		throws SyntaxException
		{
		draw (tree, new Group());
		dwg.add (group);
		}

	/**
	 * Draw a picture of the given tree string. The picture is added to the
	 * given {@linkplain edu.rit.draw.item.Group} object.
	 *
	 * @param  tree   Tree as a string in Newick Standard format.
	 * @param  group  Drawing group object.
	 *
	 * @exception  SyntaxException
	 *     Thrown if there was a syntax error in <TT>tree</TT>.
	 */
	public void draw
		(String tree,
		 Group group)
		throws SyntaxException
		{
		// Initialize.
		this.treestring = tree;
		this.tipCount = 0;
		this.group = group;

		// Set up matcher.
		matcher = pattern.matcher (treestring);
		nextToken();

		// Parse the tree, building up the drawing as it goes.
		Point root = parseNode();
		if (! ";".equals (token)) syntaxError ("; expected");
		nextToken();
		if (token != null) syntaxError ("Extra characters");

		// Add final line to root.
		group.append (new Line().to(root).hby(-H));
		}

// Hidden operations.

	/**
	 * Draw a node.
	 *
	 * @param  tree   DNA sequence tree.
	 * @param  index  Node index.
	 *
	 * @return  Point to which to attach line from parent.
	 */
	private Point drawNode
		(DnaSequenceTree tree,
		 int index)
		{
		Point p;

		// Draw a tip node.
		boolean isTip = tree.child1 (index) == -1;
		if (isTip)
			{
			p = new Point (0, tipCount * V);
			++ tipCount;
			}

		// Draw an interior node.
		else
			{
			Point p1 = drawNode (tree, tree.child1 (index));
			Point p2 = drawNode (tree, tree.child2 (index));
			double x = Math.min (p1.x(), p2.x()) - H;
			double y1 = p1.y();
			double y2 = p2.y();
			p = new Point (x, (y1 + y2)/2);
			group.append (new Line().to(p1).hto(x).vto(p2).hto(p2));
			}

		// Draw node name.
		DnaSequence seq = tree.seq (index);
		if (seq != null) drawName (seq.name(), p);

		// Draw branch length.
		Double brlen = tree.branchLength (index);
		if (brlen != null) drawBranchLength (brlen, p, isTip);

		return p;
		}

	/**
	 * Draw a node name.
	 *
	 * @param  name  Node name, or null if none.
	 * @param  p     Node's anchor point.
	 */
	private void drawName
		(String name,
		 Point p)
		{
		if (name != null)
			{
			group.append (new Text().text(name).w(p.e(GAP)));
			}
		}

	/**
	 * Draw a branch length.
	 *
	 * @param  brlen  Branch length.
	 * @param  p      Node's anchor point.
	 * @param  isTip  True if node is a tip node, false otherwise.
	 */
	private void drawBranchLength
		(double brlen,
		 Point p,
		 boolean isTip)
		{
		if (isTip)
			{
			group.append
				(new Text().text(format.format(brlen)).se(p.n(VGAP).w(GAP)));
			}
		else
			{
			group.append
				(new Text().text(format.format(brlen)).s(p.n(VGAP).w(H/2)));
			}
		}

	/**
	 * Parse a node.
	 *
	 * @return  Point to which to attach line from parent.
	 *
	 * @exception  SyntaxException
	 *     Thrown if a syntax error occurred.
	 */
	private Point parseNode()
		throws SyntaxException
		{
		if (! "(".equals (token)) syntaxError ("( expected");
		nextToken();
		Point p = parseChildList();
		if (! ")".equals (token)) syntaxError (") expected");
		nextToken();
		parseNameLength (p, false);
		return p;
		}

	/**
	 * Parse a name and/or length.
	 *
	 * @param  p      Point at which to anchor name and branch length.
	 * @param  isTip  True if node is a tip node, false otherwise.
	 *
	 * @exception  SyntaxException
	 *     Thrown if a syntax error occurred.
	 */
	private void parseNameLength
		(Point p,
		 boolean isTip)
		throws SyntaxException
		{
		if (isName (token))
			{
			drawName (token.replaceAll ("_", " "), p);
			nextToken();
			}
		if (":".equals (token))
			{
			nextToken();
			double brlen = 0.0;
			try
				{
				brlen = new Double (token);
				}
			catch (NullPointerException exc)
				{
				syntaxError ("Branch length expected");
				}
			catch (NumberFormatException exc)
				{
				syntaxError ("Branch length expected");
				}
			nextToken();
			drawBranchLength (brlen, p, isTip);
			}
		}

	/**
	 * Determine if the given token is a name.
	 *
	 * @param  token  Token.
	 *
	 * @return  True if <TT>token</TT> is a name, false otherwise.
	 */
	private boolean isName
		(String token)
		{
		return
			token != null &&
			! "(".equals (token) &&
			! ")".equals (token) &&
			! ",".equals (token) &&
			! ":".equals (token) &&
			! ";".equals (token);
		}

	/**
	 * Parse a child list.
	 *
	 * @return  Point to which to attach line from parent.
	 *
	 * @exception  SyntaxException
	 *     Thrown if a syntax error occurred.
	 */
	private Point parseChildList()
		throws SyntaxException
		{
		LinkedList<Point> childList = new LinkedList<Point>();
		double x = 0.0;
		double ymin = 0.0;
		double ymax = 0.0;
		for (;;)
			{
			Point p = parseChild();
			childList.add (p);
			x = Math.min (x, p.x());
			if (! ",".equals (token)) break;
			nextToken();
			}
		x -= H;
		for (Point p : childList)
			{
			group.append (new Line().to(p).hto(x));
			}
		ymin = childList.get(0).y();
		ymax = childList.get(childList.size()-1).y();
		group.append (new Line().to(x,ymin).to(x,ymax));
		return new Point (x, (ymin+ymax)/2);
		}

	/**
	 * Parse a child.
	 *
	 * @return  Point to which to attach line from parent.
	 *
	 * @exception  SyntaxException
	 *     Thrown if a syntax error occurred.
	 */
	private Point parseChild()
		throws SyntaxException
		{
		if ("(".equals (token))
			{
			return parseNode();
			}
		else
			{
			return parseTip();
			}
		}

	/**
	 * Parse a tip.
	 *
	 * @return  Point to which to attach line from parent.
	 *
	 * @exception  SyntaxException
	 *     Thrown if a syntax error occurred.
	 */
	private Point parseTip()
		throws SyntaxException
		{
		Point p = new Point (0, tipCount * V);
		++ tipCount;
		parseNameLength (p, true);
		return p;
		}

	/**
	 * Get the next token being parsed.
	 */
	private void nextToken()
		{
		if (matcher.find())
			{
			token = treestring.substring (matcher.start(), matcher.end());
			index = matcher.start();
			}
		else
			{
			token = null;
			index = treestring.length();
			}
		}

	/**
	 * Report a syntax error.
	 *
	 * @param  msg  Error message.
	 *
	 * @exception  SyntaxException
	 *     Thrown if a syntax error occurred.
	 */
	private void syntaxError
		(String msg)
		throws SyntaxException
		{
		throw new SyntaxException (treestring, index, msg);
		}

// Helper classes.

	/**
	 * Class TreeDrawing.SyntaxException is an exception thrown if there was a
	 * syntax error in a tree string.
	 *
	 * @author  Alan Kaminsky
	 * @version 15-Jul-2008
	 */
	public static class SyntaxException
		extends Exception
		{
		private String treestring;
		private int index;

		/**
		 * Construct a new syntax exception.
		 */
		private SyntaxException
			(String treestring,
			 int index,
			 String msg)
			{
			super (msg);
			this.treestring = treestring;
			this.index = index;
			}

		/**
		 * Print a syntax error message on the given print stream.
		 *
		 * @param  out  Print stream.
		 */
		public void printSyntaxError
			(PrintStream out)
			{
			synchronized (out)
				{
				out.println ("Syntax error: " + getMessage());
				out.println (treestring);
				for (int i = 0; i < index; ++ i) out.print (' ');
				out.println ('^');
				}
			}
		}

	}
