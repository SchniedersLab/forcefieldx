//******************************************************************************
//
// File:    ProteinLocalAlignmentSeq.java
// Package: edu.rit.compbio.seq
// Unit:    Class edu.rit.compbio.seq.ProteinLocalAlignmentSeq
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

/**
 * Class ProteinLocalAlignmentSeq provides an object that does local alignments
 * of {@linkplain ProteinSequence}s. For further information, see the base class
 * {@linkplain ProteinLocalAlignment}.
 * <P>
 * The <TT>align()</TT> method is designed to be executed by a single thread.
 * Thus, this class is suitable for use in a sequential program or one process
 * of a cluster parallel program.
 *
 * @author  Alan Kaminsky
 * @version 01-Jul-2008
 */
public class ProteinLocalAlignmentSeq
	extends ProteinLocalAlignment
	{

// Exported constructors.

	/**
	 * Construct a new protein sequence local alignment object.
	 */
	public ProteinLocalAlignmentSeq()
		{
		super();
		}

// Exported operations.

	/**
	 * Align the query sequence and the subject sequence.
	 *
	 * @return  Alignment.
	 */
	public Alignment align()
		{
		// Verify preconditions.
		if (A == null)
			{
			throw new IllegalStateException
				("ProteinLocalAlignmentSeq.align(): Query sequence not set");
			}
		if (B == null)
			{
			throw new IllegalStateException
				("ProteinLocalAlignmentSeq.align(): Subject sequence not set");
			}
		int M = A.length - 1;
		int N = B.length - 1;

		// Do the Smith-Waterman algorithm in a single thread.
		int maxScore = 0;
		int theQueryFinish = 0;
		int theSubjectFinish = 0;
		for (int i = 1; i <= M; ++ i)
			{
			int A_i = A[i];
			int[] delta_A_i = delta[A_i];
			int[] S_im1 = S[i-1];
			int[] S_i = S[i];
			int[] GA_im1 = GA[i-1];
			int[] GA_i = GA[i];
			int[] GB_i = GB[i];
			int B_j, S_i_j, GA_i_j, GB_i_j;
			for (int j = 1; j <= N; ++ j)
				{
				B_j = B[j];
				GA_i_j = S_im1[j] + g;
				GA_i_j = Math.max (GA_i_j, GA_im1[j] + h);
				GB_i_j = S_i[j-1] + g;
				GB_i_j = Math.max (GB_i_j, GB_i[j-1] + h);
				S_i_j = S_im1[j-1] + delta_A_i[B_j];
				S_i_j = Math.max (S_i_j, GA_i_j);
				S_i_j = Math.max (S_i_j, GB_i_j);
				S_i_j = Math.max (S_i_j, 0);
				if (S_i_j > maxScore)
					{
					maxScore = S_i_j;
					theQueryFinish = i;
					theSubjectFinish = j;
					}
				S_i[j] = S_i_j;
				GA_i[j] = GA_i_j;
				GB_i[j] = GB_i_j;
				}
			}

		// Do the traceback.
		return computeTraceback (maxScore, theQueryFinish, theSubjectFinish);
		}

	}
