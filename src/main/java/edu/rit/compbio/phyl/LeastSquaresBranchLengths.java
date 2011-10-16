//******************************************************************************
//
// File:    LeastSquaresBranchLengths.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.LeastSquaresBranchLengths
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

import edu.rit.numeric.NonNegativeLeastSquares;

/**
 * Class LeastSquaresBranchLengths provides methods for computing least squares
 * branch lengths in a {@linkplain DnaSequenceTree}.
 *
 * @author  Alan Kaminsky
 * @version 23-Jul-2008
 */
public class LeastSquaresBranchLengths
	{

// Prevent construction.

	private LeastSquaresBranchLengths()
		{
		}

// Exported operations.

	/**
	 * Compute the squared error in the given DNA sequence tree's branch
	 * lengths. When <TT>squaredError()</TT> is called:
	 * <UL>
	 * <LI>
	 * Every tip node in the <TT>tree</TT> must be associated with a DNA
	 * sequence of the same length.
	 * <LI>
	 * Every node in the <TT>tree</TT> either is associated with a branch
	 * length, or is assumed to have a branch length of 0 if not associated with
	 * a branch length.
	 * </UL>
	 * <P>
	 * This method computes the distance between each pair of DNA sequences
	 * using the given <TT>dcalc</TT> object (the "direct distance"); computes
	 * the distance between each pair of DNA sequences by adding up the branch
	 * lengths along the path through the tree between the sequences (the "tree
	 * distance"); and returns the sum of the squares of the differences between
	 * the direct distance and the tree distance for each pair of DNA sequences.
	 *
	 * @param  tree   DNA sequence tree.
	 * @param  dcalc  Object to calculate distances between DNA sequences.
	 *
	 * @return  Squared error.
	 */
	public static double squaredError
		(DnaSequenceTree tree,
		 Distance dcalc)
		{
		// Get tree information.
		double[] brlen = getBranchLengths (tree);
		int[] tip = getTipNodes (tree);
		boolean[][] rootPath = getRootPaths (tree, tip);
		int N = tip.length;

		// Scan all pairs of DNA sequences and compute squared error.
		double sqrerr = 0.0;
		for (int i = 0; i < N-1; ++ i)
			{
			DnaSequence seq_i = tree.seq (tip[i]);
			for (int j = i+1; j < N; ++ j)
				{
				DnaSequence seq_j = tree.seq (tip[j]);
				double d_direct = dcalc.distance (seq_i, seq_j);
				double d_tree = treeDistance (i, j, brlen, rootPath);
				double err = d_direct - d_tree;
				sqrerr += err*err;
				}
			}

		return sqrerr;
		}

	/**
	 * Find the least squares branch lengths for the given DNA sequence tree.
	 * When <TT>squaredError()</TT> is called, every tip node in the
	 * <TT>tree</TT> must be associated with a DNA sequence of the same length.
	 * <P>
	 * This method calculates the branch lengths such that the squared error, as
	 * defined in the <TT>squaredError()</TT> method, is minimized. Each node of
	 * the tree is associated with the calculated branch length. The squared
	 * error is returned.
	 * <P>
	 * This method uses a nonnegative linear least squares solver (class
	 * {@linkplain edu.rit.numeric.NonNegativeLeastSquares
	 * edu.rit.numeric.NonNegativeLeastSquares}) to calculate the branch
	 * lengths. Thus, all branch lengths will be nonnegative; some may be 0.
	 *
	 * @param  tree   DNA sequence tree.
	 * @param  dcalc  Object to calculate distances between DNA sequences.
	 *
	 * @return  Squared error.
	 */
	public static double solve
		(DnaSequenceTree tree,
		 Distance dcalc)
		{
		// Get tree information.
		double[] brlen = getBranchLengths (tree);
		int[] tip = getTipNodes (tree);
		boolean[][] rootPath = getRootPaths (tree, tip);
		int L = brlen.length;
		int N = tip.length;
		int P = N*(N - 1)/2;

		// Set up nonnegative least squares solver. Number of rows = number of
		// pairs of DNA sequences (P). Number of columns = number of branch
		// lengths (L).
		NonNegativeLeastSquares solver = new NonNegativeLeastSquares (P, L);

		// Scan all pairs of DNA sequences. For every pair of sequences p:
		// - Input vector b[p] = distance between the pair of sequences as
		//   returned by dcalc.distance().
		// - For every branch k:
		//   - Input matrix a[p][k] = 1 if the branch is on the path between the
		//     pair of sequences, = 0 otherwise.
		int p = 0;
		double[] a_p;
		for (int i = 0; i < N-1; ++ i)
			{
			DnaSequence seq_i = tree.seq (tip[i]);
			boolean[] rootPath_i = rootPath[i];
			for (int j = i+1; j < N; ++ j)
				{
				DnaSequence seq_j = tree.seq (tip[j]);
				boolean[] rootPath_j = rootPath[j];
				solver.b[p] = dcalc.distance (seq_i, seq_j);
				a_p = solver.a[p];
				for (int k = 0; k < L; ++ k)
					{
					a_p[k] = rootPath_i[k] ^ rootPath_j[k] ? 1.0 : 0.0;
					}
				++ p;
				}
			}

		// Find the solution.
		solver.solve();

		// Store branch lengths back in tree (except root has no branch length).
		for (int i = 0; i < L; ++ i)
			{
			tree.branchLength (i, solver.x[i]);
			}
		tree.branchLength (tree.root(), null);

		// Return squared error.
		return solver.normsqr;
		}

// Hidden operations.

	/**
	 * Get the branch lengths from the given tree.
	 *
	 * @param  tree  DNA sequence tree.
	 *
	 * @return  Array of branch lengths, indexed by tree node index.
	 */
	private static double[] getBranchLengths
		(DnaSequenceTree tree)
		{
		int L = tree.length();
		double[] brlen = new double [L];
		for (int i = 0; i < L; ++ i)
			{
			Double b = tree.branchLength(i);
			if (b != null) brlen[i] = b;
			}
		return brlen;
		}

	/**
	 * Get the indexes of the tip nodes in the given tree.
	 *
	 * @param  tree  DNA sequence tree.
	 *
	 * @return  Array of tip node indexes.
	 */
	private static int[] getTipNodes
		(DnaSequenceTree tree)
		{
		int L = tree.length();
		int N = (L + 1)/2;
		int[] tip = new int [N];
		int j = 0;
		for (int i = 0; i < L; ++ i)
			{
			if (tree.child1(i) == -1)
				{
				tip[j++] = i;
				}
			}
		return tip;
		}

	/**
	 * Get the paths from the tip nodes to the root node in the given tree.
	 *
	 * @param  tree  DNA sequence tree.
	 * @param  tip   Array of tip node indexes.
	 *
	 * @return  Array of paths. First index is tip node number. Second index is
	 *          branch number. The element at indexes <TT>[t,b]</TT> is true if
	 *          branch <TT>b</TT> is on the path from tip number <TT>t</TT> to
	 *          the root, false otherwise.
	 */
	private static boolean[][] getRootPaths
		(DnaSequenceTree tree,
		 int[] tip)
		{
		int L = tree.length();
		int N = tip.length;
		boolean[][] rootPath = new boolean [N] [L];
		for (int i = 0; i < N; ++ i)
			{
			boolean[] rootPath_i = rootPath[i];
			int j = tip[i];
			while (j != -1)
				{
				rootPath_i[j] = true;
				j = tree.parent(j);
				}
			}
		return rootPath;
		}

	/**
	 * Compute the tree distance between the given tip nodes.
	 *
	 * @param  tip1      First tip node number.
	 * @param  tip2      Second tip node number.
	 * @param  brlen     Array of branch lengths.
	 * @param  rootPath  Array of paths from tips to root.
	 *
	 * @return  Tree distance.
	 */
	private static double treeDistance
		(int tip1,
		 int tip2,
		 double[] brlen,
		 boolean[][] rootPath)
		{
		boolean[] rootPath1 = rootPath[tip1];
		boolean[] rootPath2 = rootPath[tip2];

		// For each branch that is the path from tip1 to root or the path from
		// tip2 to root but not both, add up the branch lengths.
		int L = rootPath[0].length;
		double d = 0.0;
		for (int i = 0; i < L; ++ i)
			{
			if (rootPath1[i] ^ rootPath2[i]) d += brlen[i];
			}

		return d;
		}

	}
