//*****************************************************************************
//
// File:    MaximumParsimonyResults.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.MaximumParsimonyResults
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

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

import java.util.ConcurrentModificationException;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Class MaximumParsimonyResults contains the results of a maximum parsimony
 * phylogenetic tree construction algorithm. The results include:
 * <UL>
 * <LI>
 * List of maximum parsimony phylogenetic trees. The list contains zero or more
 * tree signatures representing the phylogenetic trees found during the search.
 * If no tree with a score better than the initial bound was found, the list is
 * empty. For further information about the tree signatures, see the
 * <TT>toTree()</TT> method of class {@linkplain DnaSequenceList}.
 * <P><LI>
 * Parsimony score. This contains the best (smallest) Fitch parsimony score
 * found during the search. This is also the Fitch parsimony score of the
 * phylogenetic trees in the tree list (if any).
 * </UL>
 * <P>
 * <I>Note:</I> Class MaximumParsimonyResults is not multiple thread safe. It is
 * intended to be used as a per-thread variable in a parallel program.
 *
 * @author  Alan Kaminsky
 * @version 25-Oct-2008
 */
public class MaximumParsimonyResults
	implements Iterable<int[]>, Externalizable
	{

// Hidden constants.

	// Extra padding in tree list.
	private static final int PAD = 32;

// Hidden data members.

	// List of tree signatures, its size, and its capacity.
	private int[][] treeList;
	private int size;
	private int capacity;

	// Parsimony score.
	private int score;

	// For detecting modifications during an iteration.
	private int modCount;

	// Extra padding to avert cache interference.
	private static long p0, p1, p2, p3, p4, p5, p6, p7;
	private static long p8, p9, pa, pb, pc, pd, pe, pf;

// Exported constructors.

	/**
	 * Construct a new uninitialized maximum parsimony results object. This
	 * constructor is for use only by object deserialization.
	 */
	public MaximumParsimonyResults()
		{
		}

	/**
	 * Construct a new maximum parsimony results object. The tree list is
	 * initialized to an empty list with the given capacity. The parsimony score
	 * is initialized to <TT>Integer.MAX_VALUE</TT>.
	 *
	 * @param  capacity  Capacity.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>capacity</TT> &le; 0.
	 */
	public MaximumParsimonyResults
		(int capacity)
		{
		if (capacity <= 0)
			{
			throw new IllegalArgumentException
				("MaximumParsimonyResults(): capacity (= "+capacity+
				 ") illegal");
			}
		this.size = 0;
		this.capacity = capacity;
		this.treeList = new int [capacity+PAD] [];
		this.score = Integer.MAX_VALUE;
		}

	/**
	 * Construct a new maximum parsimony results object that is a copy of the
	 * given maximum parsimony results object.
	 *
	 * @param  results  Maximum parsimony results object to copy.
	 */
	public MaximumParsimonyResults
		(MaximumParsimonyResults results)
		{
		this.size = results.size;
		this.capacity = results.capacity;
		this.treeList = new int [results.capacity+PAD] [];
		for (int i = 0; i < results.size; ++ i)
			{
			this.treeList[i] = (int[]) results.treeList[i].clone();
			}
		this.score = results.score;
		}

// Exported operations.

	/**
	 * Clear this maximum parsimony results object. Afterwards, the tree list is
	 * empty and the parsimony score is <TT>Integer.MAX_VALUE</TT>.
	 */
	public void clear()
		{
		++ modCount;
		for (int i = 0; i < size; ++ i) treeList[i] = null;
		size = 0;
		score = Integer.MAX_VALUE;
		}

	/**
	 * Add the given tree with the given parsimony score to this maximum
	 * parsimony results object. The following invariant is maintained: This
	 * maximum parsimony results object contains only those trees with the
	 * smallest parsimony score seen so far; and only the first <I>C</I> such
	 * trees are stored, where <I>C</I> is the capacity.
	 *
	 * @param  tree   Tree signature.
	 * @param  score  Tree's parsimony score.
	 */
	public void add
		(int[] tree,
		 int score)
		{
		++ modCount;
		if (score < this.score)
			{
			clear();
			this.score = score;
			}
		if (score == this.score && size < capacity)
			{
			treeList[size] = (int[]) tree.clone();
			++ size;
			}
		}

	/**
	 * Add all the trees in the given maximum parsimony results object to this
	 * maximum parsimony results object. The following invariant is maintained:
	 * This maximum parsimony results object contains only those trees with the
	 * smallest parsimony score seen so far; and only the first <I>C</I> such
	 * trees are stored, where <I>C</I> is the capacity.
	 *
	 * @param  results  Maximum parsimony results object containing trees to
	 *                  add.
	 */
	public void addAll
		(MaximumParsimonyResults results)
		{
		++ modCount;
		if (results.score < this.score)
			{
			clear();
			this.score = results.score;
			}
		if (results.score == this.score)
			{
			int i = 0;
			while (this.size < this.capacity && i < results.size)
				{
				this.treeList[this.size] = (int[]) results.treeList[i].clone();
				++ this.size;
				++ i;
				}
			}
		}

	/**
	 * Returns the size of this maximum parsimony results object.
	 *
	 * @return  Size (number of trees stored).
	 */
	public int size()
		{
		return size;
		}

	/**
	 * Returns the capacity of this maximum parsimony results object.
	 *
	 * @return  Capacity (maximum number of trees that can be stored).
	 */
	public int capacity()
		{
		return capacity;
		}

	/**
	 * Returns the tree at the given index in this maximum parsimony results
	 * object.
	 *
	 * @param  i  Index, 0 &le; <TT>i</TT> &le; <TT>size()-1</TT>.
	 *
	 * @return  Tree signature at index <TT>i</TT>.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public int[] tree
		(int i)
		{
		if (0 > i || i >= size)
			{
			throw new IndexOutOfBoundsException
				("MaximumParsimonyScore.tree(): i (= "+i+") out of bounds");
			}
		return treeList[i];
		}

	/**
	 * Returns the parsimony score of this maximum parsimony results object.
	 *
	 * @return  Parsimony score.
	 */
	public int score()
		{
		return this.score;
		}

	/**
	 * Sets the parsimony score of this maximum parsimony results object.
	 *
	 * @param  score  Parsimony score.
	 */
	public void score
		(int score)
		{
		this.score = score;
		}

	/**
	 * Reduce this maximum parsimony results object's score to the given score.
	 * If this object's score is less than or equal to <TT>score</TT>, this
	 * object is unchanged. If this object's score is greater than
	 * <TT>score</TT>, this object is cleared and its score is set to
	 * <TT>score</TT>.
	 *
	 * @param  score  Parsimony score.
	 */
	public void reduceScore
		(int score)
		{
		if (this.score > score)
			{
			clear();
			this.score = score;
			}
		}

	/**
	 * Get an iterator for the trees in this maximum parsimony results object.
	 * The iterator does not support removing elements. The iterator is a
	 * "fail-fast" iterator that throws a ConcurrentModificationException if one
	 * thread changes this maximum parsimony results object while another thread
	 * is iterating over the trees.
	 *
	 * @return  Iterator.
	 */
	public Iterator<int[]> iterator()
		{
		return new Iterator<int[]>()
			{
			private int originalModCount = modCount;
			private int i;
			public boolean hasNext()
				{
				if (modCount != originalModCount)
					{
					throw new ConcurrentModificationException();
					}
				return i < size;
				}
			public int[] next()
				{
				if (modCount != originalModCount)
					{
					throw new ConcurrentModificationException();
					}
				if (i >= size)
					{
					throw new NoSuchElementException();
					}
				return treeList[i++];
				}
			public void remove()
				{
				throw new UnsupportedOperationException();
				}
			};
		}

	/**
	 * Write this maximum parsimony results object to the given object output
	 * stream.
	 *
	 * @param  out  Object output stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void writeExternal
		(ObjectOutput out)
		throws IOException
		{
		out.writeInt (size);
		out.writeInt (capacity);
		out.writeInt (score);
		for (int i = 0; i < size; ++ i)
			{
			int[] treeList_i = treeList[i];
			int n = treeList_i.length;
			out.writeInt (n);
			for (int j = 0; j < n; ++ j)
				{
				out.writeInt (treeList_i[j]);
				}
			}
		}

	/**
	 * Read this maximum parsimony results object from the given object input
	 * stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException
		{
		size = in.readInt();
		capacity = in.readInt();
		score = in.readInt();
		treeList = new int [capacity+PAD] [];
		for (int i = 0; i < size; ++ i)
			{
			int n = in.readInt();
			int[] treeList_i = new int [n];
			treeList[i] = treeList_i;
			for (int j = 0; j < n; ++ j)
				{
				treeList_i[j] = in.readInt();
				}
			}
		}

	}
