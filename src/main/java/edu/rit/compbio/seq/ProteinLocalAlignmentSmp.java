//******************************************************************************
//
// File:    ProteinLocalAlignmentSmp.java
// Package: edu.rit.compbio.seq
// Unit:    Class edu.rit.compbio.seq.ProteinLocalAlignmentSmp
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

import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import edu.rit.util.Range;

/**
 * Class ProteinLocalAlignmentSmp provides an object that does local alignments
 * of {@linkplain ProteinSequence}s. For further information, see the base class
 * {@linkplain ProteinLocalAlignment}.
 * <P>
 * The <TT>align()</TT> method is designed to be executed by a {@linkplain
 * edu.rit.pj.ParallelTeam} of threads. Thus, this class is suitable for use in
 * an SMP parallel program or a hybrid parallel program.
 * <P>
 * As an example of how the computation is performed in parallel while obeying
 * the sequential dependencies in the Smith-Waterman algorithm, suppose the
 * query sequence has 100 elements, the subject sequence has 500 elements, and
 * the parallel team has 4 threads. The 500 columns of the scoring matrix
 * <I>S</I> are partitioned equally among the threads: thread 0 gets columns
 * 1..125, thread 1 gets columns 126..250, thread 2 gets columns 251..375,
 * thread 3 gets columns 376..500. Then <I>S</I> is computed in parallel in a
 * series of rounds:
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top"><I>Round:</I>&nbsp;&nbsp;</TD>
 * <TD ALIGN="left" VALIGN="top"><I>Thread 0 computes:</I>&nbsp;&nbsp;</TD>
 * <TD ALIGN="left" VALIGN="top"><I>Thread 1 computes:</I>&nbsp;&nbsp;</TD>
 * <TD ALIGN="left" VALIGN="top"><I>Thread 2 computes:</I>&nbsp;&nbsp;</TD>
 * <TD ALIGN="left" VALIGN="top"><I>Thread 3 computes:</I>&nbsp;&nbsp;</TD>
 * </TR>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">1</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[1][1..125]</TD>
 * <TD ALIGN="left" VALIGN="top">&mdash;</TD>
 * <TD ALIGN="left" VALIGN="top">&mdash;</TD>
 * <TD ALIGN="left" VALIGN="top">&mdash;</TD>
 * </TR>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">2</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[2][1..125]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[1][126..250]</TD>
 * <TD ALIGN="left" VALIGN="top">&mdash;</TD>
 * <TD ALIGN="left" VALIGN="top">&mdash;</TD>
 * </TR>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">3</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[3][1..125]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[2][126..250]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[1][251..375]</TD>
 * <TD ALIGN="left" VALIGN="top">&mdash;</TD>
 * </TR>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">4</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[4][1..125]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[3][126..250]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[2][251..375]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[1][376..500]</TD>
 * </TR>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">5</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[5][1..125]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[4][126..250]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[3][251..375]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[2][376..500]</TD>
 * </TR>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">&hellip;</TD>
 * <TD ALIGN="left" VALIGN="top">&hellip;</TD>
 * <TD ALIGN="left" VALIGN="top">&hellip;</TD>
 * <TD ALIGN="left" VALIGN="top">&hellip;</TD>
 * <TD ALIGN="left" VALIGN="top">&hellip;</TD>
 * </TR>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">99</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[99][1..125]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[98][126..250]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[97][251..375]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[96][376..500]</TD>
 * </TR>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">100</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[100][1..125]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[99][126..250]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[98][251..375]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[97][376..500]</TD>
 * </TR>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">101</TD>
 * <TD ALIGN="left" VALIGN="top">&mdash;</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[100][126..250]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[99][251..375]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[98][376..500]</TD>
 * </TR>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">102</TD>
 * <TD ALIGN="left" VALIGN="top">&mdash;</TD>
 * <TD ALIGN="left" VALIGN="top">&mdash;</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[100][251..375]</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[99][376..500]</TD>
 * </TR>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">103</TD>
 * <TD ALIGN="left" VALIGN="top">&mdash;</TD>
 * <TD ALIGN="left" VALIGN="top">&mdash;</TD>
 * <TD ALIGN="left" VALIGN="top">&mdash;</TD>
 * <TD ALIGN="left" VALIGN="top"><I>S</I>[100][376..500]</TD>
 * </TR>
 * </TABLE>
 * <P>
 * After a short startup period, all columns of <I>S</I> are being computed in
 * parallel, with different threads working on different rows so as to obey the
 * sequential dependencies. For example, <I>S</I>[4][126] is computed (by thread
 * 1 in round 5) after <I>S</I>[3][125] (by thread 0 in round 3),
 * <I>S</I>[3][126] (by thread 1 in round 4), and <I>S</I>[4][125] (by thread 0
 * in round 4).
 *
 * @author  Alan Kaminsky
 * @version 02-Jul-2008
 */
public class ProteinLocalAlignmentSmp
	extends ProteinLocalAlignment
	{

// Hidden data members.

	private ParallelTeam team;

// Exported constructors.

	/**
	 * Construct a new protein sequence local alignment object.
	 *
	 * @param  team  Parallel thread team that will compute the alignment.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>team</TT> is null.
	 */
	public ProteinLocalAlignmentSmp
		(ParallelTeam team)
		{
		super();
		if (team == null)
			{
			throw new NullPointerException
				("ProteinLocalAlignmentSmp(): team is null");
			}
		this.team = team;
		}

// Exported operations.

	/**
	 * Align the query sequence and the subject sequence. The parallel thread
	 * team specified to the constructor computes the alignment in parallel.
	 *
	 * @return  Alignment.
	 *
	 * @exception  Exception
	 *     Thrown if an error occurred.
	 */
	public Alignment align()
		throws Exception
		{
		// Verify preconditions.
		if (A == null)
			{
			throw new IllegalStateException
				("ProteinLocalAlignmentSmp.align(): Query sequence not set");
			}
		if (B == null)
			{
			throw new IllegalStateException
				("ProteinLocalAlignmentSmp.align(): Subject sequence not set");
			}
		final int M = A.length - 1;
		final int N = B.length - 1;
		final int K = team.getThreadCount();
		final int lastRound = M + K - 1;

		// Initialize global finish point reduction variable.
		final FinishPoint gblfp = new FinishPoint();
		gblfp.maxScore = 0;
		gblfp.theQueryFinish = 0;
		gblfp.theSubjectFinish = 0;

		// Do the Smith-Waterman algorithm in the parallel thread team.
		team.execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				int threadIndex = getThreadIndex();

				// Determine range of columns for this thread.
				Range range = new Range (1, N) .subrange (K, threadIndex);
				int jlb = range.lb();
				int jub = range.ub();

				// Initialize per-thread finish point.
				int maxScore = 0;
				int theQueryFinish = 0;
				int theSubjectFinish = 0;

				// Do all rounds.
				for (int round = 1; round <= lastRound; ++ round)
					{
					// Row for this thread in this round is round number offset
					// by thread index. If row is out of bounds, do nothing this
					// round.
					int i = round - threadIndex;
					if (1 <= i && i <= M)
						{
						int A_i = A[i];
						int[] delta_A_i = delta[A_i];
						int[] S_im1 = S[i-1];
						int[] S_i = S[i];
						int[] GA_im1 = GA[i-1];
						int[] GA_i = GA[i];
						int[] GB_i = GB[i];
						int B_j, S_i_j, GA_i_j, GB_i_j;

						// Do only this thread's columns.
						for (int j = jlb; j <= jub; ++ j)
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

					// Wait for all threads to complete this round.
					barrier();
					}

				// After all rounds, reduce per-thread finish point into global
				// finish point.
				gblfp.setToBest (maxScore, theQueryFinish, theSubjectFinish);
				}
			});

		// Do the traceback in a single thread, starting from global finish
		// point.
		return computeTraceback
			(gblfp.maxScore, gblfp.theQueryFinish, gblfp.theSubjectFinish);
		}

// Hidden helper classes.

	// A record of information about the local alignment finish point.
	private static class FinishPoint
		{
		// Alignment score.
		public int maxScore;

		// Query sequence index.
		public int theQueryFinish;

		// Subject sequence index.
		public int theSubjectFinish;

		// Set this finish point to the best of itself and the given finish
		// point. Multiple thread safe method.
		public synchronized void setToBest
			(int maxScore,
			 int theQueryFinish,
			 int theSubjectFinish)
			{
			if ((maxScore > this.maxScore) ||
				(maxScore == this.maxScore &&
					theQueryFinish < this.theQueryFinish) ||
				(maxScore == this.maxScore &&
					theQueryFinish == this.theQueryFinish &&
					theSubjectFinish < this.theSubjectFinish))
				{
				this.maxScore = maxScore;
				this.theQueryFinish = theQueryFinish;
				this.theSubjectFinish = theSubjectFinish;
				}
			}
		}

	}
