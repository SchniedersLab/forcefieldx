//******************************************************************************
//
// File:    ProteinLocalAlignment.java
// Package: edu.rit.compbio.seq
// Unit:    Class edu.rit.compbio.seq.ProteinLocalAlignment
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

import java.io.ByteArrayOutputStream;

/**
 * Class ProteinLocalAlignment is the abstract base class for an object that
 * does local alignments of {@linkplain ProteinSequence}s. Class
 * ProteinLocalAlignment uses the Smith-Waterman algorithm with a substitution
 * matrix (such as BLOSUM-62) and affine gap penalties. To do a local alignment:
 * <OL TYPE=1>
 * <LI>
 * Create an instance of class ProteinLocalAlignment.
 * <P><LI>
 * If desired, call <TT>setSubstitutionMatrix()</TT> to set the substitution
 * matrix. If not set, the default is the BLOSUM-62 substitution matrix.
 * <P><LI>
 * If desired, call <TT>setGapExistencePenalty()</TT> to set the gap existence
 * penalty. If not set, the default is &minus;11.
 * <P><LI>
 * If desired, call <TT>setGapExtensionPenalty()</TT> to set the gap extension
 * penalty. If not set, the default is &minus;1.
 * <P><LI>
 * Call <TT>setQuerySequence()</TT> to set the query sequence.
 * <P><LI>
 * Call <TT>setSubjectSequence()</TT> to set the subject sequence.
 * <P><LI>
 * Call <TT>align()</TT> to do the alignment between the query sequence and the
 * subject sequence in a single thread. The <TT>align()</TT> method returns an
 * {@linkplain Alignment} object.
 * </OL>
 * <P>
 * To do another local alignment with the same query sequence and a different
 * subject sequence, repeat steps 6-7. To do another local alignment with a
 * different query sequence, repeat steps 5-7.
 * <P>
 * Subclasses of class ProteinLocalAlignment implement the Smith-Waterman
 * algorithm differently.
 *
 * @author  Alan Kaminsky
 * @version 07-Jul-2008
 */
public abstract class ProteinLocalAlignment
	{

// Hidden data members.

	// Substitution matrix.
	int[][] delta = Blosum62.matrix;

	// Gap existence and extension penalties.
	int g = -11;
	int h = -1;

	// Query sequence, ID, length.
	byte[] A;
	long myQueryId;
	int myQueryLength;

	// Subject sequence, ID, length.
	byte[] B;
	long mySubjectId;
	int mySubjectLength;

	// Score matrix.
	int[][] S;

	// Gap score matrices.
	int[][] GA;
	int[][] GB;

	// Extra padding to avert cache interference.
	long p0, p1, p2, p3, p4, p5, p6, p7;
	long p8, p9, pa, pb, pc, pd, pe, pf;

// Exported constructors.

	/**
	 * Construct a new protein sequence local alignment object.
	 */
	public ProteinLocalAlignment()
		{
		}

// Exported operations.

	/**
	 * Set the protein substitution matrix. If not set, the default is the
	 * BLOSUM-62 substitution matrix.
	 * <P>
	 * The <TT>matrix</TT> must be a 27&times;27-element matrix of integers. The
	 * expression <TT>matrix[x][y]</TT> is the score when amino acid <TT>x</TT>
	 * is aligned with amino acid <TT>y</TT>, where <TT>x</TT> and <TT>y</TT>
	 * are amino acids in the range 0..27 from a {@linkplain ProteinSequence}.
	 *
	 * @param  matrix  Protein substition matrix.
	 */
	public void setSubstitutionMatrix
		(int[][] matrix)
		{
		this.delta = matrix;
		}

	/**
	 * Set the gap existence penalty. If not set, the default is &minus;11.
	 * <P>
	 * The gap existence penalty is added to the alignment score for the first
	 * position of a gap.
	 *
	 * @param  g  Gap existence penalty. Assumed to be a negative integer.
	 */
	public void setGapExistencePenalty
		(int g)
		{
		this.g = g;
		}

	/**
	 * Set the gap extension penalty. If not set, the default is &minus;1.
	 * <P>
	 * The gap extension penalty is added to the alignment score for the second
	 * and subsequent positions of a gap.
	 *
	 * @param  h  Gap extension penalty. Assumed to be a negative integer.
	 */
	public void setGapExtensionPenalty
		(int h)
		{
		this.h = h;
		}

	/**
	 * Set the query sequence. The query sequence ID is an arbitrary integer;
	 * for example, it could be the index of a sequence in a {@linkplain
	 * ProteinDatabase}.
	 *
	 * @param  theSequence  Query sequence.
	 * @param  theId        Query sequence ID.
	 */
	public void setQuerySequence
		(ProteinSequence theSequence,
		 long theId)
		{
		A = theSequence.sequence();
		myQueryId = theId;
		myQueryLength = theSequence.length();
		int M = A.length;
		if (S == null || S.length < M+32) // Extra padding
			{
			S = new int [M+32] [];
			GA = new int [M+32] [];
			GB = new int [M+32] [];
			}
		}

	/**
	 * Set the subject sequence. The subject sequence ID is an arbitrary
	 * integer; for example, it could be the index of a sequence in a
	 * {@linkplain ProteinDatabase}.
	 *
	 * @param  theSequence  Subject sequence.
	 * @param  theId        Subject sequence ID.
	 */
	public void setSubjectSequence
		(ProteinSequence theSequence,
		 long theId)
		{
		if (A == null)
			{
			throw new IllegalStateException
				("ProteinLocalAlignment.setSubjectSequence(): Query sequence not set");
			}
		B = theSequence.sequence();
		mySubjectId = theId;
		mySubjectLength = theSequence.length();
		int N = B.length;
		if (S[0] == null || S[0].length < N+32) // Extra padding
			{
			int M = S.length-32;
			for (int i = 0; i < M; ++ i)
				{
				S[i] = new int [N+32];
				GA[i] = new int [N+32];
				GB[i] = new int [N+32];
				}
			}
		}

	/**
	 * Align the query sequence and the subject sequence.
	 * <P>
	 * Calling the returned {@linkplain Alignment} object's
	 * <TT>getQueryId()</TT> method will return the query sequence ID supplied
	 * to the <TT>setQuerySequence()</TT> method. Calling the returned
	 * {@linkplain Alignment} object's <TT>getSubjectId()</TT> method will
	 * return the subject sequence ID supplied to the
	 * <TT>setSubjectSequence()</TT> method.
	 *
	 * @return  Alignment.
	 *
	 * @exception  Exception
	 *     Thrown if an error occurred.
	 */
	public abstract Alignment align()
		throws Exception;

//	/**
//	 * Dump a section of the S, GA, GB, and delta matrices to the standard
//	 * output. This method is intended for debugging.
//	 *
//	 * @param  qlb  Query lower bound index.
//	 * @param  qub  Query upper bound index.
//	 * @param  slb  Subject lower bound index.
//	 * @param  sub  Subject upper bound index.
//	 */
//	public void dump
//		(int qlb,
//		 int qub,
//		 int slb,
//		 int sub)
//		{
//		System.out.println ("S");
//		for (int i = qlb; i <= qub; ++ i)
//			{
//			int[] S_i = S[i];
//			for (int j = slb; j <= sub; ++ j)
//				{
//				if (j > slb) System.out.print ('\t');
//				System.out.print (S_i[j]);
//				}
//			System.out.println();
//			}
//		System.out.println ("GA");
//		for (int i = qlb; i <= qub; ++ i)
//			{
//			int[] GA_i = GA[i];
//			for (int j = slb; j <= sub; ++ j)
//				{
//				if (j > slb) System.out.print ('\t');
//				System.out.print (GA_i[j]);
//				}
//			System.out.println();
//			}
//		System.out.println ("GB");
//		for (int i = qlb; i <= qub; ++ i)
//			{
//			int[] GB_i = GB[i];
//			for (int j = slb; j <= sub; ++ j)
//				{
//				if (j > slb) System.out.print ('\t');
//				System.out.print (GB_i[j]);
//				}
//			System.out.println();
//			}
//		System.out.println ("Delta");
//		for (int i = qlb; i <= qub; ++ i)
//			{
//			int A_i = A[i];
//			for (int j = slb; j <= sub; ++ j)
//				{
//				if (j > slb) System.out.print ('\t');
//				System.out.print (delta[A_i][B[j]]);
//				}
//			System.out.println();
//			}
//		}

// Hidden operations.

	/**
	 * Compute the traceback and return the resulting alignment.
	 *
	 * @param  theScore          Alignment score.
	 * @param  theQueryFinish    Location of alignment finish point in query
	 *                           sequence.
	 * @param  theSubjectFinish  Location of alignment finish point in subject
	 *                           sequence.
	 *
	 * @return  Alignment.
	 */
	Alignment computeTraceback
		(int theScore,
		 int theQueryFinish,
		 int theSubjectFinish)
		{
		// Set up alignment object.
		Alignment alignment = new Alignment();
		alignment.myQueryId = this.myQueryId;
		alignment.mySubjectId = this.mySubjectId;
		alignment.myQueryLength = this.myQueryLength;
		alignment.mySubjectLength = this.mySubjectLength;

		// Special case: No alignment found.
		if (theScore == 0)
			{
			alignment.myTraceback = new byte [0];
			return alignment;
			}

		// For recording alignment state at each position.
		ByteArrayOutputStream baos = new ByteArrayOutputStream();

		// Trace backwards until we reach a score of 0.
		int i = theQueryFinish;
		int j = theSubjectFinish;
		int theQueryStart = i;
		int theSubjectStart = j;
		int state = 0;
		while (S[i][j] != 0)
			{
			switch (state)
				{
				case 0: // Tracing back through table S
					if (S[i][j] == S[i-1][j-1] + delta[A[i]][B[j]])
						{
						baos.write (Alignment.QUERY_ALIGNED_WITH_SUBJECT);
						theQueryStart = i;
						theSubjectStart = j;
						-- i;
						-- j;
						}
					else if (S[i][j] == GA[i][j])
						{
						state = 1;
						}
					else
						{
						state = 2;
						}
					break;
				case 1: // Tracing back through table GA
					baos.write (Alignment.QUERY_ALIGNED_WITH_GAP);
					if (GA[i][j] == S[i-1][j] + g) state = 0;
					theQueryStart = i;
					-- i;
					break;
				case 2: // Tracing back through table GB
					baos.write (Alignment.SUBJECT_ALIGNED_WITH_GAP);
					if (GB[i][j] == S[i][j-1] + g) state = 0;
					theSubjectStart = j;
					-- j;
					break;
				}
			}

		// Record results.
		alignment.myScore = theScore;
		alignment.myQueryStart = theQueryStart;
		alignment.mySubjectStart = theSubjectStart;
		alignment.myQueryFinish = theQueryFinish;
		alignment.mySubjectFinish = theSubjectFinish;
		alignment.myTraceback = baos.toByteArray();
		return alignment;
		}

	}
