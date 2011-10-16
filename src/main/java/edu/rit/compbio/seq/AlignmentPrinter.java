//******************************************************************************
//
// File:    AlignmentPrinter.java
// Package: edu.rit.compbio.seq
// Unit:    Class edu.rit.compbio.seq.AlignmentPrinter
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

package edu.rit.compbio.seq;

import java.io.PrintStream;

import java.util.Formatter;

/**
 * Class AlignmentPrinter provides an object that prints an {@linkplain
 * Alignment}.
 *
 * @author  Alan Kaminsky
 * @version 01-Jul-2008
 */
public class AlignmentPrinter
	{

// Hidden data members.

	private PrintStream out;
	private AlignmentStats stats;
	private int[][] delta = Blosum62.matrix;

// Exported constructors.

	/**
	 * Construct a new alignment printer.
	 *
	 * @param  out    Print stream on which to print (e.g.,
	 *                <TT>System.out</TT>).
	 * @param  stats  Object for calculating alignment statistics.
	 */
	public AlignmentPrinter
		(PrintStream out,
		 AlignmentStats stats)
		{
		this.out = out;
		this.stats = stats;
		}

// Exported operations.

	/**
	 * Set the substitution matrix. If not set, the default is the BLOSUM-62
	 * substitution matrix.
	 * <P>
	 * The expression <TT>matrix[x][y]</TT> is the score when sequence element
	 * <TT>x</TT> is aligned with sequence element <TT>y</TT>.
	 *
	 * @param  matrix  Substition matrix.
	 */
	public void setSubstitutionMatrix
		(int[][] matrix)
		{
		this.delta = matrix;
		}

	/**
	 * Print a summary of the given alignment. The printout consists of:
	 * <UL>
	 * <LI>
	 * Subject sequence description (columns 1-60).
	 * <LI>
	 * Bit score (columns 63-67).
	 * <LI>
	 * <I>E</I>-value (columns 70-75).
	 * </UL>
	 *
	 * @param  alignment  Alignment.
	 * @param  subject    Subject sequence.
	 */
	public void printSummary
		(Alignment alignment,
		 Sequence subject)
		{
		String desc = subject.description();
		if (desc.length() > 60) desc = desc.substring(0,57) + "...";
		double bitScore = stats.bitScore (alignment);
		double eValue = stats.eValue (alignment);
		out.format
			("%-60s  %5s  %s%n",
			 desc,
			 formatBitScore (bitScore),
			 formatEValue (eValue));
		}

	/**
	 * Print details of the given alignment. The printout consists of:
	 * <UL>
	 * <LI>
	 * Subject sequence description and length.
	 * <LI>
	 * Bit score, raw score, and <I>E</I>-value.
	 * <LI>
	 * Number and percent of identities, positives, and gaps.
	 * <LI>
	 * The aligned query and subject sequences.
	 * </UL>
	 *
	 * @param  alignment  Alignment.
	 * @param  query      Query sequence.
	 * @param  subject    Subject sequence.
	 */
	public void printDetails
		(Alignment alignment,
		 Sequence query,
		 Sequence subject)
		{
		// Print subject sequence description and length.
		out.println (subject.description());
		out.format ("Length = %d%n", subject.length());
		out.println();

		// Print bit score, raw score, and <I>E</I>-value.
		double bitScore = stats.bitScore (alignment);
		double rawScore = stats.rawScore (alignment);
		double eValue = stats.eValue (alignment);
		out.format
			("Score = %s bits (%.0f), Expect = %s%n",
			 formatBitScore (bitScore),
			 rawScore,
			 formatEValue (eValue));

		// Count identities, positives, and gaps.
		int identities = 0;
		int positives = 0;
		int gaps = 0;
		int qi = alignment.myQueryStart;
		int si = alignment.mySubjectStart;
		int n = alignment.myTraceback.length;
		for (int i = 0; i < n; ++ i)
			{
			switch (alignment.myTraceback[n-1-i])
				{
				case Alignment.QUERY_ALIGNED_WITH_SUBJECT:
					byte query_qi = query.mySequence[qi];
					byte subject_si = subject.mySequence[si];
					if (query_qi == subject_si)
						{
						++ identities;
						++ positives;
						}
					else if (delta[query_qi][subject_si] > 0)
						{
						++ positives;
						}
					++ qi;
					++ si;
					break;
				case Alignment.QUERY_ALIGNED_WITH_GAP:
					++ gaps;
					++ qi;
					break;
				case Alignment.SUBJECT_ALIGNED_WITH_GAP:
					++ gaps;
					++ si;
					break;
				}
			}

		// Print identities, positives, and gaps.
		out.format
			("Identities = %d/%d (%.0f%%), ",
			 identities,
			 n,
			 ((double) identities)/((double) n)*100.0);
		out.format
			("Positives = %d/%d (%.0f%%), ",
			 positives,
			 n,
			 ((double) positives)/((double) n)*100.0);
		out.format
			("Gaps = %d/%d (%.0f%%)%n",
			 gaps,
			 n,
			 ((double) gaps)/((double) n)*100.0);
		out.println();

		// Print aligned sequences.
		qi = alignment.myQueryStart;
		si = alignment.mySubjectStart;
		int i = 0;
		int j;
		int qj;
		int sj;
		while (i < n)
			{
			qj = qi;
			out.format ("Query%6d  ", qj);
			j = 0;
			while (j < 60 && i+j < n)
				{
				switch (alignment.myTraceback[n-1-i-j])
					{
					case Alignment.QUERY_ALIGNED_WITH_SUBJECT:
					case Alignment.QUERY_ALIGNED_WITH_GAP:
						out.print (query.charAt(qj));
						++ qj;
						break;
					case Alignment.SUBJECT_ALIGNED_WITH_GAP:
						out.print ('-');
						break;
					}
				++ j;
				}
			out.format ("%6d%n", qj-1);
			qj = qi;
			sj = si;
			out.print ("             ");
			j = 0;
			while (j < 60 && i+j < n)
				{
				switch (alignment.myTraceback[n-1-i-j])
					{
					case Alignment.QUERY_ALIGNED_WITH_SUBJECT:
						byte query_qj = query.mySequence[qj];
						byte subject_sj = subject.mySequence[sj];
						if (query_qj == subject_sj)
							{
							out.print (query.charAt(qj));
							}
						else if (delta[query_qj][subject_sj] > 0)
							{
							out.print ('+');
							}
						else
							{
							out.print (' ');
							}
						++ qj;
						++ sj;
						break;
					case Alignment.QUERY_ALIGNED_WITH_GAP:
						out.print (' ');
						++ qj;
						break;
					case Alignment.SUBJECT_ALIGNED_WITH_GAP:
						out.print (' ');
						++ sj;
						break;
					}
				++ j;
				}
			out.println();
			qj = qi;
			sj = si;
			out.format ("Sbjct%6d  ", sj);
			j = 0;
			while (j < 60 && i+j < n)
				{
				switch (alignment.myTraceback[n-1-i-j])
					{
					case Alignment.QUERY_ALIGNED_WITH_SUBJECT:
						out.print (subject.charAt(sj));
						++ qj;
						++ sj;
						break;
					case Alignment.QUERY_ALIGNED_WITH_GAP:
						out.print ('-');
						++ qj;
						break;
					case Alignment.SUBJECT_ALIGNED_WITH_GAP:
						out.print (subject.charAt(sj));
						++ sj;
						break;
					}
				++ j;
				}
			out.format ("%6d%n", sj-1);
			out.println();
			qi = qj;
			si = sj;
			i += j;
			}
		}

// Hidden operations.

	/**
	 * Format the given bit score.
	 *
	 * @param  bitScore  Bit score.
	 *
	 * @return  Formatted bit score.
	 */
	private static String formatBitScore
		(double bitScore)
		{
		Formatter f = new Formatter();
		if (bitScore >= 100.0) f.format ("%.0f", bitScore);
		else if (bitScore >= 10.0) f.format ("%.1f", bitScore);
		else f.format ("%.2f", bitScore);
		return f.toString();
		}

	/**
	 * Format the given <I>E</I>-value.
	 *
	 * @param  eValue  <I>E</I>-value.
	 *
	 * @return  Formatted <I>E</I>-value.
	 */
	private static String formatEValue
		(double eValue)
		{
		Formatter f = new Formatter();
		if (eValue < 1.0e-199) f.format ("0.0");
		else if (eValue < 0.001) f.format ("%.0e", eValue);
		else if (eValue < 0.1) f.format ("%.3f", eValue);
		else if (eValue < 1.0) f.format ("%.2f", eValue);
		else if (eValue < 10.0) f.format ("%.1f", eValue);
		else f.format ("%.0f", eValue);
		return f.toString();
		}

	}
