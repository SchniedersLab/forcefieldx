//*****************************************************************************
//
// File:    MaximumParsimonyExhSeq.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.MaximumParsimonyExhSeq
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

/**
 * Class MaximumParsimonyExhSeq provides a sequential algorithm for maximum
 * parsimony phylogenetic tree construction using exhaustive search.
 * <P>
 * To perform a search, the program must:
 * <OL TYPE=1>
 * <LI>
 * Create an instance of class MaximumParsimonyExhSeq, passing in a {@linkplain
 * DnaSequenceList} of the DNA sequences in the tree and a {@linkplain
 * MaximumParsimonyResults} object to hold the search results.
 * <P><LI>
 * Call the <TT>findTrees()</TT> method. The results of the search are returned
 * in the {@linkplain MaximumParsimonyResults} object specified to the
 * constructor.
 * </OL>
 * <P>
 * <I>Note:</I> Class MaximumParsimonyExhSeq is not multiple thread safe; it is
 * intended to be used in a single-threaded program.
 *
 * @author  Alan Kaminsky
 * @version 21-Nov-2008
 */
public class MaximumParsimonyExhSeq
	{

// Hidden data members.

	// List of DNA sequences with which to construct trees.
	private DnaSequenceList seqList;

	// For holding search results.
	private MaximumParsimonyResults results;

	// Length of each DNA sequence.
	private int L;

	// Number of DNA sequences.
	private int N;

	// Tree capacity.
	private int C;

	// Stack of DNA sequence trees.
	private DnaSequenceTree[] treeStack;

	// Stack of auxiliary DNA sequence arrays.
	DnaSequence[][] seqArrayStack;

	// Tree signature currently being searched.
	private int[] signature;

// Exported constructors.

	/**
	 * Construct a new maximum parsimony phylogenetic tree construction
	 * algorithm object.
	 *
	 * @param  seqList  DNA sequence list.
	 * @param  results  Object in which to store the results.
	 */
	public MaximumParsimonyExhSeq
		(DnaSequenceList seqList,
		 MaximumParsimonyResults results)
		{
		// Record parameters.
		this.seqList = seqList;
		this.results = results;

		// Initialize.
		L = seqList.seq(0).length();
		N = seqList.length();
		C = 2*N - 1;

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
		signature = new int [N];
		}

// Exported operations.

	/**
	 * Find the maximum parsimony phylogenetic tree(s) in the search graph. The
	 * DNA sequence list was specified to the constructor. The results are
	 * stored in the {@linkplain MaximumParsimonyResults} object specified to
	 * the constructor.
	 */
	public void findTrees()
		{
		// Initialize tree signature as specified by <startLevel> and <vertex1>.
		signature[0] = 0;
		for (int i = 1; i < N; ++ i)
			{
			signature[i] = -1;
			}

		// Traverse all levels of the search graph.
		int level = 1;
		boolean done = false;
		results.clear();
		while (! done)
			{
			DnaSequenceTree prevTree = treeStack[level-1];

			// If we have reached the bottom of the search graph, we have a
			// tentative solution.
			if (level == N)
				{
				int tentativeScore = prevTree.seq (prevTree.root()) .score();

				// Record tentative solution.
				results.add (signature, tentativeScore);

				// Go to previous level.
				-- level;
				if (level == 1) done = true;
				}

			// If there are no more positions to try at this level, reset
			// position at this level and go to previous level.
			else if (signature[level] == 2*(level - 1))
				{
				signature[level] = -1;
				-- level;
				if (level == 1) done = true;
				}

			// If there are more positions to try at this level, add the DNA
			// sequence to the tree at the next position and go to next level.
			else
				{
				++ signature[level];
				DnaSequenceTree currTree = treeStack[level];
				currTree.copy (prevTree);
				int tip = currTree.add (signature[level], seqList.seq(level));
				FitchParsimony.updateScore
					(currTree, tip, seqArrayStack[level]);
				++ level;
				}
			}
		}

	}
