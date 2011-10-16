//*****************************************************************************
//
// File:    Results.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.Results
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

import edu.rit.swing.DisplayableIO;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import java.util.Arrays;
import java.util.Date;
import java.util.Formatter;
import java.util.Map;

/**
 * Class Results provides a method for reporting the results of a maximum
 * parsimony phylogenetic tree construction program. The results are stored in
 * files in a specified directory. The following files are created:
 * <UL>
 * <LI>
 * <TT>"index.html"</TT> has all the results. Open this file in a web browser to
 * view the results.
 * <P><LI>
 * <TT>"tree_001.png"</TT>, <TT>"tree_002.png"</TT>, and so on have pictures of
 * the phylogenetic trees. Each tip node is labeled with a DNA sequence name.
 * Each interior node is labeled with the number of state changes at that node.
 * <P><LI>
 * <TT>"sequences.phy"</TT> has a copy of the input DNA sequences in interleaved
 * PHYLIP format.
 * <P><LI>
 * <TT>"trees.txt"</TT> has all the phylogenetic trees, one per line, in Newick
 * Standard format.
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 25-Oct-2008
 */
public class Results
	{

// Prevent construction.

	private Results()
		{
		}

// Exported operations.

	/**
	 * Report the results of a maximum parsimony phylogenetic tree construction
	 * program.
	 *
	 * @param  directory
	 *     Directory in which to store the files. If the directory does not
	 *     exist, it is created.
	 * @param  programName
	 *     Name of the program.
	 * @param  hostName
	 *     Name of the host computer, or null if none.
	 * @param  K
	 *     Number of parallel threads, or 0 if none.
	 * @param  infile
	 *     Input file.
	 * @param  originalSeqList
	 *     List of DNA sequences read from the input file, in the order in which
	 *     the sequences appeared in the input file.
	 * @param  sortedSeqList
	 *     List of DNA sequences read from the input file, in the order in which
	 *     the sequences were added to the phylogenetic trees during the search.
	 * @param  initialBound
	 *     Initial bound for branch-and-bound search, or 0 if no initial bound.
	 * @param  treeStoreLimit
	 *     Maximum number of trees to store.
	 * @param  results
	 *     Maximum parsimony results object containing results of the search.
	 * @param  t1
	 *     System clock snapshot at beginning.
	 * @param  t2
	 *     System clock snapshot after preprocessing.
	 * @param  t3
	 *     System clock snapshot after computation.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void report
		(File directory,
		 String programName,
		 String hostName,
		 int K,
		 File infile,
		 DnaSequenceList originalSeqList,
		 DnaSequenceList sortedSeqList,
		 int initialBound,
		 int treeStoreLimit,
		 MaximumParsimonyResults results,
		 long t1,
		 long t2,
		 long t3)
		throws IOException
		{
		// If directory does not exist, create it.
		if (! directory.exists()) directory.mkdirs();

		// Begin writing "index.html" file.
		File indexfile = new File (directory, "index.html");
		PrintStream out =
			new PrintStream
				(new BufferedOutputStream
					(new FileOutputStream (indexfile)));
		out.println ("<HTML>");
		out.println ("<HEAD>");
		out.print   ("<TITLE>Maximum Parsimony Phylogenetic Tree Search Results -- ");
		out.print   (infile);
		out.println ("</TITLE>");
		out.println ("<STYLE TYPE=\"text/css\">");
		out.println ("<!--");
		out.println ("* {font-family: Arial, Helvetica, Sans-Serif;}");
		out.println ("body {font-size: small;}");
		out.println ("h1 {font-size: 140%; font-weight: bold;}");
		out.println ("h2 {font-size: 120%; font-weight: bold;}");
		out.println ("h3 {font-size: 100%; font-weight: bold;}");
		out.println ("table {font-size: 100%;}");
		out.println ("tt {font-family: Courier, Monospace; font-size: 100%;}");
		out.println ("tt b {font-family: Courier, Monospace; font-size: 100%; font-weight: normal; background: #e8e8e8;}");
		out.println ("pre {font-family: Courier, Monospace; font-size: 100%;}");
		out.println ("pre b {font-family: Courier, Monospace; font-size: 100%; font-weight: normal; background: #e8e8e8;}");
		out.println ("-->");
		out.println ("</STYLE>");
		out.println ("</HEAD>");
		out.println ("<BODY>");

		// Write overview section of "index.html" file.
		out.println ("<H1>Maximum Parsimony Phylogenetic Tree Search Results</H1>");
		out.println ("<P><HR/><H2>Summary</H2></P>");
		out.println ("<P>");
		out.println ("<TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Date/time:&nbsp;&nbsp;</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (new Date (t3));
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Program:&nbsp;&nbsp;</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (programName);
		out.println ("</TD>");
		out.println ("</TR>");
		if (hostName != null)
			{
			out.println ("<TR>");
			out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Computer:&nbsp;&nbsp;</TD>");
			out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
			out.print   (hostName.replaceAll("<","&lt;").replaceAll(">","&gt;"));
			out.println ("</TD>");
			out.println ("</TR>");
			}
		if (K > 0)
			{
			out.println ("<TR>");
			out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Parallel threads:&nbsp;&nbsp;</TD>");
			out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
			out.print   (K);
			out.println ("</TD>");
			out.println ("</TR>");
			}
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Sequence file:&nbsp;&nbsp;</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (infile);
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Number of sequences:&nbsp;&nbsp;</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (originalSeqList.length());
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Number of sites:&nbsp;&nbsp;</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (originalSeqList.seq(0).length());
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Number of informative sites:&nbsp;&nbsp;</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (originalSeqList.informativeSiteCount());
		out.println ("</TD>");
		out.println ("</TR>");
		if (initialBound >= 0)
			{
			out.println ("<TR>");
			out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Parsimony score initial bound:&nbsp;&nbsp;</TD>");
			out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
			out.print   (initialBound);
			out.println ("</TD>");
			out.println ("</TR>");
			}
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Maximum trees saved:&nbsp;&nbsp;</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (treeStoreLimit);
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Number of trees:&nbsp;&nbsp;</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (results.size());
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Number of state changes:&nbsp;&nbsp;</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (results.score());
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("</TABLE>");
		out.println ("</P>");
		out.println ("<P>");
		out.println ("<A HREF=\"#sequences\">Sequences</A>");
		out.println ("<BR/><A HREF=\"#distances\">Distances</A>");
		out.println ("<BR/><A HREF=\"#trees\">Trees</A>");
		out.println ("</P>");

		// Write a copy of the input sequences to "sequences.phy" file.
		originalSeqList.write
			(new File (directory, "sequences.phy"), 70, true, false);

		// Write sequences section of "index.html" file.
		out.println ("<P><HR/><A NAME=\"sequences\"><H2>Sequences</H2></A></P>");
		out.println ("<P>Sequences in interleaved PHYLIP format:&nbsp;&nbsp;<A HREF=\"sequences.phy\">sequences.phy</A></P>");
		out.println ("<P>");
		out.println ("Notes:");
		out.print   ("<BR/>&diams;&nbsp;&nbsp;Number of sequences: ");
		out.println (originalSeqList.length());
		out.print   ("<BR/>&diams;&nbsp;&nbsp;Number of sites: ");
		out.println (originalSeqList.seq(0).length());
		out.print   ("<BR/>&diams;&nbsp;&nbsp;Number of informative sites: ");
		out.println (originalSeqList.informativeSiteCount());
		out.println ("<BR/>&diams;&nbsp;&nbsp;Informative sites <TT><B>marked</B></TT>");
		out.println ("</P>");
		out.println ("<PRE>");
		originalSeqList.write (out, 70, true, true);
		out.println ("</PRE>");

		// Write distances section of "index.html" file.
		out.println ("<P><HR/><A NAME=\"distances\"><H2>Distances</H2></A></P>");
		Distance dcalc_H = new HammingDistance();
		printDistanceMatrix
			(out, "Hamming distances:", "%.0f", originalSeqList, dcalc_H);
		Distance dcalc_JC = new JukesCantorDistance();
		printDistanceMatrix
			(out, "Jukes-Cantor distances:", "%.2f", originalSeqList, dcalc_JC);

		// Write trees section of "index.html" file.
		out.println ("<P><HR/><A NAME=\"trees\"><H2>Trees</H2></A></P>");
		out.println ("<P>Trees in Newick Standard format:&nbsp;&nbsp;<A HREF=\"trees.txt\">trees.txt</A></P>");
		out.println ("<P>");
		out.println ("Notes:");
		out.println ("<BR/>&diams;&nbsp;&nbsp;These are unrooted trees");
		out.println ("<BR/>&diams;&nbsp;&nbsp;Each interior node marked with number of state changes at that node");
		out.println ("<BR/>&diams;&nbsp;&nbsp;Each branch marked with least squares branch length, Jukes-Cantor distances");
		out.print   ("<BR/>&diams;&nbsp;&nbsp;Number of trees: ");
		out.print   (results.size());
		out.print   ("<BR/>&diams;&nbsp;&nbsp;Total number of state changes: ");
		out.print   (results.score());
		out.println ("</P>");

		// Make an array of all trees in ascending order of squared error and
		// average root height.
		int L = results.size();
		TreeInfo[] allTrees = new TreeInfo [L];
		for (int i = 0; i < L; ++ i)
			{
			// Reconstruct tree.
			DnaSequenceTree tree =
				sortedSeqList.toTree (results.tree(i));

			// Compute parsimony score at each node.
			FitchParsimony.computeScore (tree);

			// Compute least squares branch lengths.
			double sqrerr = LeastSquaresBranchLengths.solve (tree, dcalc_JC);

			// Compute average root height.
			double avgheight = averageRootHeight (tree);

			// Add to array.
			allTrees[i] = new TreeInfo (tree, sqrerr, avgheight);
			}
		Arrays.sort (allTrees);

		// Write all trees.
		PrintStream treeout =
			new PrintStream
				(new BufferedOutputStream
					(new FileOutputStream
						(new File (directory, "trees.txt"))));
		Drawing drawing = Drawing.defaultDrawing();
		TreeDrawing artist = new TreeDrawing();
		for (int i = 0; i < L; ++ i)
			{
			DnaSequenceTree tree = allTrees[i].tree;
			double sqrerr = allTrees[i].sqrerr;
			double avgheight = allTrees[i].avgheight;

			// Print tree to "index.html" file.
			out.print   ("<P><H3>Tree ");
			out.print   (i + 1);
			out.print   (" of ");
			out.print   (L);
			out.print   ("</H3></P>");
			out.println ("<PRE>");
			out.println (tree);
			out.println ("</PRE>");

			// Print tree to "trees.txt" file.
			treeout.println (tree);

			// Write tree drawing to "tree_XXX.png" file.
			drawing.clear();
			artist.draw (tree);
			StringBuilder builder = new StringBuilder();
			Formatter formatter = new Formatter (builder);
			formatter.format ("tree_%03d.png", i + 1);
			String treefilename = builder.toString();
			DisplayableIO.writeGrayscalePNGFile
				(drawing, new File (directory, treefilename));

			// Display tree drawing in "index.html" file.
			out.print   ("<IMG SRC=\"");
			out.print   (treefilename);
			out.println ("\"/>");
			out.println ("<P>");
			out.print   ("Squared error = ");
			out.println (sqrerr);
			out.println ("</P>");
			}
		treeout.close();

		// Finish writing "index.html" file.
		long t4 = System.currentTimeMillis();
		out.println ("<P><HR/></P>");
		out.println ("<P>");
		out.println ("<TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>");
		out.println ("<TR>");
		out.print   ("<TD ALIGN=\"right\" VALIGN=\"top\">");
		out.print   (t2 - t1);
		out.println ("&nbsp;</TD>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">msec preprocessing</TD>");
		out.println ("</TR>");
		out.println ("<TR>");
		out.print   ("<TD ALIGN=\"right\" VALIGN=\"top\">");
		out.print   (t3 - t2);
		out.println ("&nbsp;</TD>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">msec calculation</TD>");
		out.println ("</TR>");
		out.println ("<TR>");
		out.print   ("<TD ALIGN=\"right\" VALIGN=\"top\">");
		out.print   (t4 - t3);
		out.println ("&nbsp;</TD>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">msec postprocessing</TD>");
		out.println ("</TR>");
		out.println ("<TR>");
		out.print   ("<TD ALIGN=\"right\" VALIGN=\"top\">");
		out.print   (t4 - t1);
		out.println ("&nbsp;</TD>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">msec total</TD>");
		out.println ("</TR>");
		out.println ("</TABLE>");
		out.println ("</P>");
		out.println ("<P>");
		out.println ("<TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.println ("Powered by Parallel Java:&nbsp;&nbsp;");
		out.println ("</TD>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.println ("<A HREF=\"http://www.cs.rit.edu/~ark/pj.shtml\">http://www.cs.rit.edu/~ark/pj.shtml</A>");
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.println ("Developed by Alan Kaminsky:&nbsp;&nbsp;");
		out.println ("</TD>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.println ("<A HREF=\"http://www.cs.rit.edu/~ark/\">http://www.cs.rit.edu/~ark/</A>");
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("</TABLE>");
		out.println ("</P>");
		out.println ("</BODY>");
		out.println ("</HTML>");
		out.close();
		}

// Hidden operations.

	/**
	 * Print a distance matrix.
	 *
	 * @param  out      Where to print.
	 * @param  title    Title.
	 * @param  format   Format string for printing distances.
	 * @param  seqList  DNA sequence list.
	 * @param  dcalc    Object to calculate distances between sequences.
	 */
	private static void printDistanceMatrix
		(PrintStream out,
		 String title,
		 String format,
		 DnaSequenceList seqList,
		 Distance dcalc)
		{
		out.print   ("<P>");
		out.print   (title);
		out.println ("</P>");
		out.println ("<TABLE BORDER=1 CELLPADDING=2 CELLSPACING=0>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"center\" BGCOLOR=\"#e8e8e8\">&nbsp;</TD>");
		for (int j = 0; j < seqList.length(); ++ j)
			{
			DnaSequence seq_j = seqList.seq(j);
			out.print   ("<TD ALIGN=\"center\" BGCOLOR=\"#e8e8e8\">");
			out.print   (seq_j.name());
			out.println ("</TD>");
			}
		out.println ("</TR>");
		for (int i = 0; i < seqList.length(); ++ i)
			{
			DnaSequence seq_i = seqList.seq(i);
			out.println ("<TR>");
			out.print   ("<TD ALIGN=\"center\" BGCOLOR=\"#e8e8e8\">");
			out.print   (seq_i.name());
			out.println ("</TD>");
			for (int j = 0; j < seqList.length(); ++ j)
				{
				DnaSequence seq_j = seqList.seq(j);
				out.print   ("<TD ALIGN=\"center\">");
				out.format  (format, dcalc.distance (seq_i, seq_j));
				out.println ("</TD>");
				}
			out.println ("</TR>");
			}
		out.println ("</TABLE>");
		}

	/**
	 * Compute average of branch lengths from every tip node to the root in the
	 * given tree.
	 *
	 * @param  tree  DNA sequence tree.
	 */
	private static double averageRootHeight
		(DnaSequenceTree tree)
		{
		double sum = 0.0;
		int L = tree.length();
		int N = 0;
		for (int i = 0; i < L; ++ i)
			{
			if (tree.child1(i) == -1)
				{
				sum += rootHeight (tree, i);
				++ N;
				}
			}
		return sum/N;
		}

	/**
	 * Compute branch length from the given tip node to the root in the given
	 * tree.
	 *
	 * @param  tree  DNA sequence tree.
	 * @param  i   Tip node index.
	 */
	private static double rootHeight
		(DnaSequenceTree tree,
		 int i)
		{
		double sum = 0.0;
		while (i != -1)
			{
			Double brlen = tree.branchLength(i);
			if (brlen != null) sum += brlen;
			i = tree.parent(i);
			}
		return sum;
		}

// Hidden helper classes.

	/**
	 * A record of information about a DNA sequence tree.
	 */
	private static class TreeInfo
		implements Comparable<TreeInfo>
		{
		public DnaSequenceTree tree;
		public double sqrerr;
		public double avgheight;

		public TreeInfo
			(DnaSequenceTree tree,
			 double sqrerr,
			 double avgheight)
			{
			this.tree = tree;
			this.sqrerr = sqrerr;
			this.avgheight = avgheight;
			}

		public int compareTo
			(TreeInfo info)
			{
			if (this.sqrerr < info.sqrerr) return -1;
			else if (this.sqrerr > info.sqrerr) return +1;
			else if (this.avgheight < info.avgheight) return -1;
			else if (this.avgheight > info.avgheight) return +1;
			else return 0;
			}
		}

	}
