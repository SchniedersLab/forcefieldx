//******************************************************************************
//
// File:    DefaultAlignmentStats.java
// Package: edu.rit.compbio.seq
// Unit:    Class edu.rit.compbio.seq.DefaultAlignmentStats
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

/**
 * Class DefaultAlignmentStats provides an object that computes
 * statistics of an {@linkplain Alignment}. Methods are provided to compute the
 * raw score, the bit score, and the <I>E</I>-value. The formulas for these
 * statistics assume that the alignment:
 * <UL>
 * <LI>
 * is for two protein sequences,
 * <LI>
 * was calculated by the Smith-Waterman local alignment algorithm,
 * <LI>
 * using the BLOSUM-62 protein substitution matrix,
 * <LI>
 * using a gap existence penalty of &minus;11 and a gap extension penalty of
 * &minus;1,
 * <LI>
 * matching a query sequence against a database of subject sequences, where the
 * sum of the subject sequence lengths is supplied as a constructor parameter.
 * </UL>
 * <P>
 * The formulas for the bit score and <I>E</I>-value are:
 * <CENTER>
 * <I>S</I>' = (<I>&lambda;</I> <I>S</I> &minus; ln <I>K</I>)/(ln 2)
 * <BR>
 * <I>E</I> = <I>K</I> <I>m</I> <I>n</I> exp(&minus;<I>&lambda;</I> <I>S</I>)
 * </CENTER>
 * where <I>S</I> is the raw score, <I>m</I> is the query sequence length,
 * <I>n</I> is the total subject sequence length, and the parameters are
 * <I>K</I> = 0.035, <I>&lambda;</I> = 0.252.
 * <P>
 * These formulas were taken from:
 * <UL>
 * <LI>
 * <A HREF="http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html" TARGET="_top">http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html</A>
 * <LI>
 * <A HREF="http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-3.html" TARGET="_top">http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-3.html</A>
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 01-Jul-2008
 */
public class DefaultAlignmentStats
	implements AlignmentStats
	{

// Hidden constants.

	private static final double K = 0.035;
	private static final double lambda = 0.252;
	private static final double LN_K = Math.log(K);
	private static final double LN_2 = Math.log(2.0);

// Hidden data members.

	// Total subject sequence length.
	private long n;

// Exported constructors.

	/**
	 * Construct a new default alignment statistics object.
	 *
	 * @param  n  Sum of the lengths of the subject sequences in the database.
	 */
	public DefaultAlignmentStats
		(long n)
		{
		this.n = n;
		}

// Exported operations.

	/**
	 * Returns the raw score for the given alignment. A larger raw score
	 * signifies a greater degree of similarity between the query sequence and
	 * subject sequence that were aligned.
	 *
	 * @param  alignment  Alignment.
	 *
	 * @return  Raw score.
	 */
	public double rawScore
		(Alignment alignment)
		{
		return alignment.myScore;
		}

	/**
	 * Returns the bit score for the given alignment. A larger bit score
	 * signifies a greater degree of similarity between the query sequence and
	 * subject sequence that were aligned.
	 * <P>
	 * The bit score is the raw score normalized to units of "bits." Bit scores
	 * for different alignment procedures may be compared, whereas raw
	 * (unnormalized) scores for different alignment procedures may not be
	 * compared.
	 *
	 * @param  alignment  Alignment.
	 *
	 * @return  Bit score.
	 */
	public double bitScore
		(Alignment alignment)
		{
		double S = alignment.myScore;
		return (lambda*S - LN_K)/LN_2;
		}

	/**
	 * Returns the <I>E</I>-value (expect value) for the given alignment. A
	 * smaller <I>E</I>-value signifies a more statistically significant degree
	 * of similarity between the query sequence and subject sequence that were
	 * aligned.
	 * <P>
	 * The <I>E</I>-value is the expected number of alignments with a score
	 * greater than or equal to the <TT>alignment</TT>'s score when a
	 * randomly-chosen query of the same length as the query that produced the
	 * <TT>alignment</TT> is matched against the database.
	 *
	 * @param  alignment  Alignment.
	 *
	 * @return  Bit score.
	 */
	public double eValue
		(Alignment alignment)
		{
		int m = alignment.myQueryLength;
		double S = alignment.myScore;
		return K*m*n*Math.exp(-lambda*S);
		}

	/**
	 * Print information about this alignment statistics object on the given
	 * print stream.
	 *
	 * @param  out  Print stream.
	 */
	public void print
		(PrintStream out)
		{
		out.println ("K: "+K);
		out.println ("Lambda: "+lambda);
		out.println ("Matrix: BLOSUM-62");
		out.println ("Gap Penalties: Existence: -11, Extension: -1");
		}

	}
