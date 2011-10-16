//******************************************************************************
//
// File:    Upgma2.java
// Package: edu.rit.compbio.phyl.test
// Unit:    Class edu.rit.compbio.phyl.test.Upgma2
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

package edu.rit.compbio.phyl.test;

import edu.rit.compbio.phyl.DnaSequence;
import edu.rit.compbio.phyl.DnaSequenceTree;

/**
 * Class Upgma2 is a unit test program for the Unweighted Pair Group Method with
 * Arithmetic mean (UPGMA) algorithm for phylogenetic tree construction. It
 * constructs the tree for a specific distance matrix (the one on Felsenstein,
 * page 163). For further information, see:
 * <UL>
 * <LI>
 * R. Sokal and C. Michener. A statistical method for evaluating systematic
 * relationships. <I>University of Kansas Science Bulletin,</I> 38:1409-1438,
 * 1958.
 * <LI>
 * J. Felsenstein. <I>Inferring Phylogenies.</I> Sinauer Associates, 2004, pages
 * 161-166.
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 15-Jul-2008
 */
public class Upgma2
	{

// Prevent construction.

	private Upgma2()
		{
		}

// Fixed distance matrix.

	static final int ORIG_N = 8;

	static final double[][] ORIG_D = new double[][]
		{
		{  0,  32,  48,  51,  50,  48,  98, 148,   0},
		{ 32,   0,  26,  34,  29,  33,  84, 136,   0},
		{ 48,  26,   0,  42,  44,  44,  92, 152,   0},
		{ 51,  34,  42,   0,  44,  38,  86, 142,   0},
		{ 50,  29,  44,  44,   0,  24,  89, 142,   0},
		{ 48,  33,  44,  38,  24,   0,  90, 142,   0},
		{ 98,  84,  92,  86,  89,  90,   0, 148,   0},
		{148, 136, 152, 142, 142, 142, 148,   0,   0},
		{  0,   0,   0,   0,   0,   0,   0,   0,   0},
		};

	static final DnaSequence[] ORIG_SEQ = new DnaSequence[]
		{
		new DnaSequence (0, 0, "DOG"),
		new DnaSequence (0, 0, "BEAR"),
		new DnaSequence (0, 0, "RACCOON"),
		new DnaSequence (0, 0, "WEASEL"),
		new DnaSequence (0, 0, "SEAL"),
		new DnaSequence (0, 0, "SEA LION"),
		new DnaSequence (0, 0, "CAT"),
		new DnaSequence (0, 0, "MONKEY"),
		};

// Exported operations.

	/**
	 * Main program. Prints the tree in Newick Standard format.
	 * <P>
	 * Usage: java edu.rit.compbio.phyl.Upgma2
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		System.out.println (buildTree());
		}

// Hidden operations.

	/**
	 * Build a phylogenetic tree of the given DNA sequences.
	 *
	 * @return  Phylogenetic tree.
	 */
	private static DnaSequenceTree buildTree()
		{
		// Get initial DNA sequences and put each one in its own tree.
		int N = ORIG_N;
		DnaSequenceTree[] tree = new DnaSequenceTree [N+1];
		for (int i = 0; i < N; ++ i)
			{
			tree[i] = new DnaSequenceTree (1);
			tree[i].add (0, ORIG_SEQ[i]);
			}

		// Compute initial distance matrix.
		double[][] D = ORIG_D;

		// Set up array of group sizes n_i.
		int[] n = new int [N+1];
		for (int i = 0; i < N; ++ i)
			{
			n[i] = 1;
			}

		// Join trees until only one is left.
		while (N > 1)
			{
			// Find i and j for which D_i_j is smallest.
			double min_d = Double.POSITIVE_INFINITY;
			int min_i = 0;
			int min_j = 0;
			for (int i = 0; i < N-1; ++ i)
				{
				double[] D_i = D[i];
				for (int j = i+1; j < N; ++ j)
					{
					double d = D_i[j];
					if (d < min_d)
						{
						min_d = d;
						min_i = i;
						min_j = j;
						}
					}
				}

			// Compute node heights for trees <min_i> and <min_j>. Store as
			// branch lengths for now.
			DnaSequenceTree tree_i = tree[min_i];
			tree_i.branchLength (tree_i.root(), 0.5*min_d);
			DnaSequenceTree tree_j = tree[min_j];
			tree_j.branchLength (tree_j.root(), 0.5*min_d);

			// Join trees <min_i> and <min_j>. Add new tree to end of list.
			DnaSequenceTree newtree =
				new DnaSequenceTree (tree_i.length() + tree_j.length() + 1);
			newtree.join (tree_i, tree_j);
			tree[N] = newtree;
			int newn = n[min_i] + n[min_j];
			n[N] = newn;

			// Compute distance from new tree to every other tree.
			double w_i = ((double) n[min_i])/newn;
			double w_j = ((double) n[min_j])/newn;
			double[] D_n = D[N];
			for (int k = 0; k < N; ++ k)
				{
				double[] D_k = D[k];
				double D_n_k = w_i*D_k[min_i] + w_j*D_k[min_j];
				D_n[k] = D_n_k;
				D_k[N] = D_n_k;
				}

			// Swap row <N> with row <min_i> and swap row <N-1> with row
			// <min_j>, thus removing rows <min_i> and <min_j> from D, tree, and
			// n.
			double[] swap1 = D[min_i];
			D[min_i] = D[N];
			D[N] = swap1;
			swap1 = D[min_j];
			D[min_j] = D[N-1];
			D[N-1] = swap1;
			DnaSequenceTree swap2 = tree[min_i];
			tree[min_i] = tree[N];
			tree[N] = swap2;
			swap2 = tree[min_j];
			tree[min_j] = tree[N-1];
			tree[N-1] = swap2;
			int swap3 = n[min_i];
			n[min_i] = n[N];
			n[N] = swap3;
			swap3 = n[min_j];
			n[min_j] = n[N-1];
			n[N-1] = swap3;

			// Swap column <N> with column <min_i> and swap column <N-1> with
			// column <min_j>, thus removing columns <min_i> and <min_j> from D.
			for (int i = 0; i <= N; ++ i)
				{
				double[] D_i = D[i];
				double swap4 = D_i[min_i];
				D_i[min_i] = D_i[N];
				D_i[N] = swap4;
				swap4 = D_i[min_j];
				D_i[min_j] = D_i[N-1];
				D_i[N-1] = swap4;
				}

			// Took away two trees and added one.
			D[N] = null;
			tree[N] = null;
			-- N;
			}

		// Convert node heights to branch lengths.
		computeBranchLengths (tree[0], tree[0].root());

		// Return the last tree left standing.
		return tree[0];
		}

	/**
	 * Convert the node height stored in the given node of the given tree to a
	 * branch length.
	 *
	 * @param  tree   Tree.
	 * @param  index  Node index.
	 *
	 * @return  Original node height of the given node.
	 */
	private static double computeBranchLengths
		(DnaSequenceTree tree,
		 int index)
		{
		// Stop recursion when we reach a tip.
		if (index == -1) return 0.0;

		// Convert both children.
		double childHeight = computeBranchLengths (tree, tree.child1 (index));
		computeBranchLengths (tree, tree.child2 (index));

		// Convert this node.
		Double nodeHeight = tree.branchLength (index);
		if (nodeHeight == null)
			{
			return 0.0;
			}
		else
			{
			tree.branchLength (index, nodeHeight - childHeight);
			return nodeHeight;
			}
		}

	}
