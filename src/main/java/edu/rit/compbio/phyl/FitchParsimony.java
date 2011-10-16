//******************************************************************************
//
// File:    FitchParsimony.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.FitchParsimony
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
 * Class FitchParsimony provides the Fitch algorithm for computing the parsimony
 * score of a {@linkplain DnaSequenceTree}. For further information, see:
 * <UL>
 * <LI>
 * W. Fitch. Toward defining the course of evolution: minimum change for a
 * specified tree topology. <I>Systematic Zoology,</I> 20:406-416, 1971.
 * <LI>
 * J. Felsenstein. <I>Inferring Phylogenies.</I> Sinauer Associates, 2004, pages
 * 11-13.
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 18-Jul-2008
 */
public class FitchParsimony
	{

// Prevent construction.

	private FitchParsimony()
		{
		}

// Exported operations.

	/**
	 * Compute the Fitch parsimony score of the given DNA sequence tree. Call
	 * <TT>computeScore()</TT> to compute the score for an entire tree.
	 * <P>
	 * When <TT>computeScore()</TT> is called:
	 * <UL>
	 * <LI>
	 * Every tip node in <TT>tree</TT> must be associated with a DNA sequence
	 * whose score is 0.
	 * <LI>
	 * Every interior node in <TT>tree</TT> may or may not be associated with a
	 * DNA sequence.
	 * <LI>
	 * All DNA sequences in <TT>tree</TT> must be the same length.
	 * </UL>
	 * <P>
	 * When <TT>computeScore()</TT> returns:
	 * <UL>
	 * <LI>
	 * Every tip node's DNA sequence is unchanged.
	 * <LI>
	 * Every interior node that had not been associated with a DNA sequence is
	 * associated with a newly created DNA sequence.
	 * <LI>
	 * Every interior node's DNA sequence has been set to an intermediate value
	 * as determined by the Fitch algorithm, and the DNA sequence's name has
	 * been set to the number of state changes at that node.
	 * <LI>
	 * The root node's DNA sequence's score is the tree's Fitch parsimony score.
	 * The score is also returned.
	 * </UL>
	 *
	 * @param  tree  DNA sequence tree.
	 *
	 * @return  Fitch parsimony score of <TT>tree</TT>.
	 */
	public static int computeScore
		(DnaSequenceTree tree)
		{
		int root = tree.root();
		computeScore (tree, root);
		return tree.seq (root) .score();
		}

	/**
	 * Compute the Fitch parsimony score of the given node in the given DNA
	 * sequence tree.
	 *
	 * @param  tree   DNA sequence tree.
	 * @param  index  Node index.
	 */
	private static void computeScore
		(DnaSequenceTree tree,
		 int index)
		{
		// Stop recursion at a tip node.
		int child1 = tree.child1 (index);
		int child2 = tree.child2 (index);
		if (child1 == -1) return;

		// Compute scores of child nodes.
		computeScore (tree, child1);
		computeScore (tree, child2);

		// Associate a new DNA sequence with this node if necessary.
		DnaSequence seq1 = tree.seq (child1);
		DnaSequence seq2 = tree.seq (child2);
		DnaSequence seq = tree.seq (index);
		if (seq == null)
			{
			seq = new DnaSequence (seq1.length());
			tree.seq (index, seq);
			}

		// Set this node's DNA sequence to the Fitch ancestor of the two child
		// nodes' DNA sequences.
		seq.setFitchAncestor (seq1, seq2);
		seq.name (""+(seq.score()-seq1.score()-seq2.score()));
		}

	/**
	 * Update the Fitch parsimony score of the given DNA sequence tree, from the
	 * given tip node up to the root. Call <TT>updateScore()</TT> to re-compute
	 * the score when a tip node is added to the tree; this takes less time than
	 * re-computing the score for the entire tree.
	 * <P>
	 * When <TT>updateScore()</TT> is called:
	 * <UL>
	 * <LI>
	 * The node at index <TT>tip</TT> and its parent node must have been just
	 * added by the <TT>DnaSequenceTree.add()</TT> method.
	 * <LI>
	 * The node at index <TT>tip</TT> must be associated with a DNA sequence
	 * whose score is 0.
	 * <LI>
	 * All other nodes in <TT>tree</TT> must be associated with a DNA sequence
	 * containing a Fitch parsimony score, as computed by a previous call of
	 * <TT>computeScore()</TT> or <TT>updateScore()</TT>.
	 * <LI>
	 * All DNA sequences in <TT>tree</TT> must be the same length.
	 * </UL>
	 * <P>
	 * When <TT>computeScore()</TT> returns:
	 * <UL>
	 * <LI>
	 * Every interior node on the path from the node at index <TT>tip</TT> to
	 * the root node is associated with one of the sequences in
	 * <TT>seqarray</TT>. If there are <I>N</I> tip nodes in the tree, there
	 * must be at least <I>N</I>&minus;1 DNA sequences in <TT>seqarray</TT>.
	 * <LI>
	 * Every interior node's DNA sequence on the path from the node at index
	 * <TT>tip</TT> to the root node has been set to an intermediate value as
	 * determined by the Fitch algorithm.
	 * <LI>
	 * The root node's DNA sequence's score is the tree's Fitch parsimony score.
	 * The score is also returned.
	 * </UL>
	 *
	 * @param  tree      DNA sequence tree.
	 * @param  tip       Index of just-added tip node.
	 * @param  seqarray  Array of DNA sequences to be associated with interior
	 *                   nodes.
	 *
	 * @return  Fitch parsimony score of <TT>tree</TT>.
	 */
	public static int updateScore
		(DnaSequenceTree tree,
		 int tip,
		 DnaSequence[] seqarray)
		{
		int i = 0;

		// Update all nodes from tip's parent through root.
		DnaSequence seq = tree.seq (tip);
		int index = tree.parent (tip);
		while (index != -1)
			{
			int child1 = tree.child1 (index);
			int child2 = tree.child2 (index);
			DnaSequence seq1 = tree.seq (child1);
			DnaSequence seq2 = tree.seq (child2);
			seq = seqarray[i++];
			seq.setFitchAncestor (seq1, seq2);
			tree.seq (index, seq);
			index = tree.parent (index);
			}

		return seq.score();
		}

	}
