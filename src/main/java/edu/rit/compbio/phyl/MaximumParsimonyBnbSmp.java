//*****************************************************************************
//
// File:    MaximumParsimonyBnbSmp.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.MaximumParsimonyBnbSmp
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

import edu.rit.pj.reduction.IntegerOp;
import edu.rit.pj.reduction.SharedInteger;

/**
 * Class MaximumParsimonyBnbSmp provides an SMP parallel algorithm for maximum
 * parsimony phylogenetic tree construction using branch-and-bound search.
 * <P>
 * Class MaximumParsimonyBnbSmp is designed to be used in an SMP parallel
 * program that runs in one process with multiple threads. Each thread has its
 * own MaximumParsimonyBnbSmp instance. Each thread uses its own
 * MaximumParsimonyBnbSmp instance to search different sections of the search
 * graph concurrently.
 * <P>
 * The process has a shared variable, <I>bound</I>, that holds the best
 * parsimony score found so far (i.e., the bound for branch-and-bound search).
 * The <I>bound</I> variable is an instance of class {@linkplain
 * edu.rit.pj.reduction.SharedInteger}. All the MaximumParsimonyBnbSmp instances
 * within the process share the same <I>bound</I> variable. Whenever one thread
 * finds a phylogenetic tree with a better parsimony score, it notifies all the
 * threads by updating the <I>bound</I> variable.
 * <P>
 * To perform a search, the process must:
 * <OL TYPE=1>
 * <LI>
 * Call the static <TT>createBoundVariable()</TT> method to create the shared
 * <I>bound</I> variable.
 * </OL>
 * Then each thread in the process must:
 * <OL TYPE=1 START=2>
 * <LI>
 * Create its own instance of class MaximumParsimonyBnbSmp, passing in a
 * {@linkplain DnaSequenceList} of the DNA sequences in the tree, a reference to
 * the process's shared <I>bound</I> variable, and a {@linkplain
 * MaximumParsimonyResults} object to hold the search results.
 * <P><LI>
 * Call the <TT>findTrees()</TT> method one or more times, each time indicating
 * a section of the search graph to search. The results of searching that
 * section are accumulated into the {@linkplain MaximumParsimonyResults} object
 * specified to the constructor.
 * </OL>
 * <P>
 * <I>Note:</I> Class MaximumParsimonyBnbSmp is not multiple thread safe; it is
 * intended to be used as a per-thread variable.
 *
 * @author  Alan Kaminsky
 * @version 21-Nov-2008
 */
public class MaximumParsimonyBnbSmp
	{

// Exported static operations.

	/**
	 * Create a new shared <I>bound</I> variable.
	 *
	 * @param  initialBound  Initial bound for branch-and-bound search.
	 */
	public static SharedInteger createBoundVariable
		(int initialBound)
		{
		return new SharedInteger (initialBound);
		}

// Hidden data members.

	// List of DNA sequences with which to construct trees.
	private DnaSequenceList seqList;

	// Shared bound variable.
	private SharedInteger bound;

	// For holding search results.
	private MaximumParsimonyResults results;

	// Length of each DNA sequence.
	private int L;

	// Number of DNA sequences.
	private int N;

	// Tree capacity.
	private int C;

	// Number of absent states as each DNA sequence is added.
	private int[] absentStates;

	// Stack of DNA sequence trees.
	private DnaSequenceTree[] treeStack;

	// Stack of auxiliary DNA sequence arrays.
	DnaSequence[][] seqArrayStack;

	// Tree signature currently being searched.
	private int[] signature;

	// Extra padding to avert cache interference.
	private long p0, p1, p2, p3, p4, p5, p6, p7;
	private long p8, p9, pa, pb, pc, pd, pe, pf;

// Exported constructors.

	/**
	 * Construct a new maximum parsimony phylogenetic tree construction
	 * algorithm object.
	 *
	 * @param  seqList  DNA sequence list.
	 * @param  bound    Shared <I>bound</I> variable.
	 * @param  results  Object in which to store the results.
	 */
	public MaximumParsimonyBnbSmp
		(DnaSequenceList seqList,
		 SharedInteger bound,
		 MaximumParsimonyResults results)
		{
		// Record parameters.
		this.seqList = seqList;
		this.bound = bound;
		this.results = results;

		// Initialize.
		L = seqList.seq(0).length();
		N = seqList.length();
		C = 2*N - 1;

		// Compute number of absent states as each DNA sequence is added.
		absentStates = seqList.countAbsentStates();

		// Set up stack of DNA sequence trees.
		treeStack = new DnaSequenceTree [N];
		for (int i = 0; i < N; ++ i)
			{
			treeStack[i] = new DnaSequenceTree (C);
			}

		// Initialize DNA sequence tree at first level of the search graph.
		treeStack[0].add (0, seqList.seq(0));

		// Set up stack of auxiliary DNA sequence arrays.
		seqArrayStack = new DnaSequence [N] [];
		for (int i = 0; i < N; ++ i)
			{
			DnaSequence[] seqArray = new DnaSequence [i];
			seqArrayStack[i] = seqArray;
			for (int j = 0; j < i; ++ j)
				{
				seqArray[j] = new DnaSequence (L);
				}
			}

		// Set up tree signature.
		signature = new int [N+32]; // Extra padding
		}

// Exported operations.

	/**
	 * Find the maximum parsimony phylogenetic tree(s) in the given section of
	 * the search graph. The DNA sequence list was specified to the constructor.
	 * <P>
	 * The search will commence at level <I>L</I> of the search graph, 0 &le;
	 * <I>L</I> &le; <I>N</I>&minus;1, where <I>N</I> is the number of sequences
	 * in the DNA sequence list. Of the (2<I>L</I>&minus;1)!! vertices at level
	 * <I>L</I>, the search will start at the <I>V</I><SUB>1</SUB>-th such
	 * vertex and end at the <I>V</I><SUB>2</SUB>-th such vertex, 0 &le;
	 * <I>V</I><SUB>1</SUB> &le; <I>V</I><SUB>2</SUB> &le;
	 * (2<I>L</I>&minus;1)!!&nbsp;&minus;&nbsp;1. All search graph vertices at
	 * and below vertices <I>V</I><SUB>1</SUB> through <I>V</I><SUB>2</SUB>
	 * inclusive will be searched.
	 * <P>
	 * The results are accumulated into the {@linkplain MaximumParsimonyResults}
	 * object specified to the constructor. The <TT>findTrees()</TT> method will
	 * only find trees whose parsimony scores are less than or equal to the
	 * value of the <TT>bound</TT> variable specified to the constructor or the
	 * best bound found thereafter, whichever is smaller.
	 *
	 * @param  startLevel
	 *     <I>L</I>, the level of the search graph at which to commence the
	 *     search.
	 * @param  vertex1
	 *     <I>V</I><SUB>1</SUB>, the search graph vertex at level <I>L</I> at
	 *     which to start the search.
	 * @param  vertex2
	 *     <I>V</I><SUB>2</SUB>, the search graph vertex at level <I>L</I> at
	 *     which to end the search.
	 */
	public void findTrees
		(int startLevel,
		 int vertex1,
		 int vertex2)
		{
		// Initialize tree signature as specified by <startLevel> and <vertex1>.
		signature[0] = 0;
		int q = vertex1;
		for (int i = startLevel; i > 0; -- i)
			{
			int d = 2*i - 1;
			signature[i] = q % d - 1;
			q = q / d;
			}
		for (int i = startLevel + 1; i < N; ++ i)
			{
			signature[i] = -1;
			}

		// Traverse remaining levels of the search graph.
		int level = 1;
		boolean done = false;
		results.reduceScore (bound.get());
		while (! done)
			{
			DnaSequenceTree prevTree = treeStack[level-1];

			// If we have reached the bottom of the search graph, we have a
			// tentative solution.
			if (level == N)
				{
				int tentativeScore = prevTree.seq (prevTree.root()) .score();

				// Update best solution's score to reflect tentative solution's
				// score.
				int updatedScore =
					bound.reduce (tentativeScore, IntegerOp.MINIMUM);

				// If best solution's score is better than that of previous
				// solutions, discard previous solutions.
				results.reduceScore (updatedScore);

				// Record tentative solution.
				results.add (signature, tentativeScore);

				// Go to previous level.
				-- level;
				if (level == startLevel)
					{
					++ vertex1;
					if (vertex1 > vertex2) done = true;
					}
				}

			// If there are no more positions to try at this level, reset
			// position at this level and go to previous level.
			else if (signature[level] == 2*(level - 1))
				{
				signature[level] = -1;
				-- level;
				if (level == startLevel)
					{
					++ vertex1;
					if (vertex1 > vertex2) done = true;
					}
				}

			// If there are more positions to try at this level, add the DNA
			// sequence to the tree at the next position and do
			// branch-and-bound.
			else
				{
				++ signature[level];
				DnaSequenceTree currTree = treeStack[level];
				currTree.copy (prevTree);
				int tip = currTree.add (signature[level], seqList.seq(level));
				int partialScore =
					FitchParsimony.updateScore
						(currTree, tip, seqArrayStack[level]);

				// If partial parsimony score plus number of absent states in
				// the remaining levels is less than or equal to the best
				// solution's score, go to the next level (branch), otherwise
				// try the next choice at this level (bound).
				if
					((level <= startLevel) ||
					 (partialScore + absentStates[level] <= bound.get()))
					{
					++ level;
					}
				}
			}
		}

	}
