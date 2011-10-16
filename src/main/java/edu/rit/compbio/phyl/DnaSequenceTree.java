//******************************************************************************
//
// File:    DnaSequenceTree.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.DnaSequenceTree
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

import java.util.ArrayList;
import java.util.Collections;

/**
 * Class DnaSequenceTree encapsulates a rooted bifurcating tree of DNA
 * sequences. Each node in the tree is designated by an index from 0 to
 * <I>N</I>&minus;1, where <I>N</I> is the tree's <I>length.</I> The tree's
 * <I>capacity</I> <I>C</I> is the maximum number of nodes (specified when the
 * tree was constructed). For a tree to hold <I>M</I> tip nodes, the tree's
 * capacity must be <I>C</I> &ge; 2<I>M</I>&nbsp;&minus;&nbsp;1.
 *
 * @author  Alan Kaminsky
 * @version 15-Jul-2008
 */
public class DnaSequenceTree
	{

// Hidden constants.

	// Amount of extra padding in arrays.
	private static final int PAD = 32;

// Hidden helper classes.

	// Information about one node in the tree.
	private static class Node
		implements Comparable<Node>
		{
		// Parent node index, or -1 for the root node.
		public int parent;

		// First child node index, or -1 for a tip node.
		public int child1;

		// Second child node index, or -1 for a tip node.
		public int child2;

		// DNA sequence associated with this node, or null if none.
		public DnaSequence seq;

		// Branch length between this node and its parent, or null if none.
		public Double brlen;

		// Extra padding to avert cache interference.
		private long p0, p1, p2, p3, p4, p5, p6, p7;
		private long p8, p9, pa, pb, pc, pd, pe, pf;

		// Construct a new, empty node.
		public Node()
			{
			clear();
			}

		// Clear this node.
		public void clear()
			{
			this.parent = -1;
			this.child1 = -1;
			this.child2 = -1;
			this.seq = null;
			this.brlen = null;
			}

		// Set this node to be a copy of the given node.
		public void copy
			(Node node)
			{
			this.parent = node.parent;
			this.child1 = node.child1;
			this.child2 = node.child2;
			this.seq = node.seq;
			this.brlen = node.brlen;
			}

		// Compare this node with the given node. The ordering is descending
		// order of branch length. A node with no branch length uses a default
		// branch length of 0.
		public int compareTo
			(Node node)
			{
			double brlen1 = this.brlen == null ? 0.0 : this.brlen;
			double brlen2 = node.brlen == null ? 0.0 : node.brlen;
			if (brlen1 > brlen2) return -1;
			else if (brlen1 < brlen2) return +1;
			else return 0;
			}
		}

// Hidden data members.

	// Array of nodes.
	private Node[] myNode;

	// Tree capacity (maximum number of nodes).
	private int myCapacity;

	// Tree length (actual number of nodes).
	private int myLength;

	// Index of root node, or -1 if tree is empty.
	private int myRoot;

	// Extra padding to avert cache interference.
	private long p0, p1, p2, p3, p4, p5, p6, p7;
	private long p8, p9, pa, pb, pc, pd, pe, pf;

// Exported constructors.

	/**
	 * Construct a new DNA sequence tree with the given capacity. The new tree
	 * is initially empty.
	 *
	 * @param  C  Capacity.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>C</TT> &lt; 0.
	 */
	public DnaSequenceTree
		(int C)
		{
		if (C < 0)
			{
			throw new IllegalArgumentException
				("DnaSequenceTree(): C (= "+C+") < 0, illegal");
			}
		myNode = new Node [C + PAD];
		for (int i = 0; i < C; ++ i)
			{
			myNode[i] = new Node();
			}
		myCapacity = C;
		myLength = 0;
		myRoot = -1;
		}

// Exported operations.

	/**
	 * Returns the capacity of this tree.
	 *
	 * @return  Capacity <I>C</I> (maximum number of nodes).
	 */
	public int capacity()
		{
		return myCapacity;
		}

	/**
	 * Returns the length of this tree.
	 *
	 * @return  Length <I>N</I> (number of nodes).
	 */
	public int length()
		{
		return myLength;
		}

	/**
	 * Returns the root of this tree.
	 *
	 * @return  Index of the root node, or &minus;1 if this tree is empty.
	 */
	public int root()
		{
		return myRoot;
		}

	/**
	 * Returns the parent of the given node in this tree.
	 *
	 * @param  i  Node index, 0 &le; <TT>i</TT> &le; <I>N</I>&minus;1.
	 *
	 * @return  Index of the parent of node <TT>i</TT>, or &minus;1 if node
	 *          <TT>i</TT> is the root node.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public int parent
		(int i)
		{
		if (0 > i || i >= myLength)
			{
			throw new IndexOutOfBoundsException
				("DnaSequenceTree.parent(): Index (= "+i+") out of bounds");
			}
		return myNode[i].parent;
		}

	/**
	 * Returns the first child of the given node in this tree.
	 *
	 * @param  i  Node index, 0 &le; <TT>i</TT> &le; <I>N</I>&minus;1.
	 *
	 * @return  Index of the first child of node <TT>i</TT>, or &minus;1 if node
	 *          <TT>i</TT> is a tip node.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public int child1
		(int i)
		{
		if (0 > i || i >= myLength)
			{
			throw new IndexOutOfBoundsException
				("DnaSequenceTree.child1(): Index (= "+i+") out of bounds");
			}
		return myNode[i].child1;
		}

	/**
	 * Returns the second child of the given node in this tree.
	 *
	 * @param  i  Node index, 0 &le; <TT>i</TT> &le; <I>N</I>&minus;1.
	 *
	 * @return  Index of the second child of node <TT>i</TT>, or &minus;1 if
	 *          node <TT>i</TT> is a tip node.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public int child2
		(int i)
		{
		if (0 > i || i >= myLength)
			{
			throw new IndexOutOfBoundsException
				("DnaSequenceTree.child2(): Index (= "+i+") out of bounds");
			}
		return myNode[i].child2;
		}

	/**
	 * Returns the DNA sequence associated with the given node in this tree.
	 *
	 * @param  i  Node index, 0 &le; <TT>i</TT> &le; <I>N</I>&minus;1.
	 *
	 * @return  DNA sequence associated with node <TT>i</TT>, or null if no DNA
	 *          sequence is associated.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public DnaSequence seq
		(int i)
		{
		if (0 > i || i >= myLength)
			{
			throw new IndexOutOfBoundsException
				("DnaSequenceTree.seq(): Index (= "+i+") out of bounds");
			}
		return myNode[i].seq;
		}

	/**
	 * Set the DNA sequence associated with the given node in this tree.
	 * <P>
	 * <I>Note:</I> The tree contains a reference to (not a copy of)
	 * <TT>seq</TT>.
	 *
	 * @param  i    Node index, 0 &le; <TT>i</TT> &le; <I>N</I>&minus;1.
	 * @param  seq  DNA sequence associated with node <TT>i</TT>, or null if no
	 *              DNA sequence is associated.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public void seq
		(int i,
		 DnaSequence seq)
		{
		if (0 > i || i >= myLength)
			{
			throw new IndexOutOfBoundsException
				("DnaSequenceTree.seq(): Index (= "+i+") out of bounds");
			}
		myNode[i].seq = seq;
		}

	/**
	 * Returns the branch length associated with the given node in this tree.
	 *
	 * @param  i  Node index, 0 &le; <TT>i</TT> &le; <I>N</I>&minus;1.
	 *
	 * @return  Branch length between node <TT>i</TT> and its parent, or null if
	 *          no branch length is associated.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public Double branchLength
		(int i)
		{
		if (0 > i || i >= myLength)
			{
			throw new IndexOutOfBoundsException
				("DnaSequenceTree.branchLength(): Index (= "+i+
				 ") out of bounds");
			}
		return myNode[i].brlen;
		}

	/**
	 * Set the branch length associated with the given node in this tree.
	 *
	 * @param  i      Node index, 0 &le; <TT>i</TT> &le; <I>N</I>&minus;1.
	 * @param  brlen  Branch length between node <TT>i</TT> and its parent, or
	 *                null if no branch length is associated.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public void branchLength
		(int i,
		 Double brlen)
		{
		if (0 > i || i >= myLength)
			{
			throw new IndexOutOfBoundsException
				("DnaSequenceTree.branchLength(): Index (= "+i+
				 ") out of bounds");
			}
		myNode[i].brlen = brlen;
		}

	/**
	 * Clear this DNA sequence tree.
	 */
	public void clear()
		{
		int C = myCapacity;
		for (int i = 0; i < C; ++ i)
			{
			myNode[i].clear();
			}
		myLength = 0;
		myRoot = -1;
		}

	/**
	 * Set this DNA sequence tree to be a copy of the given tree. This tree's
	 * capacity is unchanged and must be greater than or equal to the given
	 * tree's length. This tree's length, DNA sequences, and branch lengths
	 * become the same as <TT>tree</TT>.
	 * <P>
	 * <I>Note:</I> This tree contains references to (not copies of) the DNA
	 * sequences in <TT>tree</TT>.
	 *
	 * @param  tree  DNA sequence tree.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>tree</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if this tree's capacity is less than
	 *     <TT>tree</TT>'s length.
	 */
	public void copy
		(DnaSequenceTree tree)
		{
		// Verify preconditions.
		if (this.myCapacity < tree.myLength)
			{
			throw new IllegalArgumentException
				("DnaSequenceTree.copy(): Capacity (= "+this.myCapacity+
				 ") too small");
			}
		int C = this.myCapacity;
		int N = tree.myLength;

		// Copy <tree>'s nodes to this tree.
		for (int i = 0; i < N; ++ i)
			{
			this.myNode[i].copy (tree.myNode[i]);
			}

		// Clear any unused nodes in this tree.
		for (int i = N; i < C; ++ i)
			{
			this.myNode[i].clear();
			}

		// Update this tree's length and root.
		this.myLength = N;
		this.myRoot = tree.myRoot;
		}

	/**
	 * Set this DNA sequence tree to be the join of the two given trees. This
	 * tree's capacity is unchanged and must be greater than or equal to
	 * <I>N</I><SUB>1</SUB>&nbsp;+&nbsp;<I>N</I><SUB>2</SUB>&nbsp;+&nbsp;1,
	 * where <I>N</I><SUB>1</SUB> and <I>N</I><SUB>2</SUB> are the two given
	 * trees' lengths. This tree's root node's first child becomes
	 * <TT>tree1</TT>. This tree's root node's second child becomes
	 * <TT>tree2</TT>. This tree's root node has no associated DNA sequence.
	 * <P>
	 * <I>Note:</I> This tree contains references to (not copies of) the DNA
	 * sequences in <TT>tree1</TT> and <TT>tree2</TT>.
	 * <P>
	 * <I>Note:</I> Both <TT>tree1</TT> and <TT>tree2</TT> must be different
	 * objects from this tree.
	 * <P>
	 * <I>Note:</I> This method may alter the index of this tree's root node.
	 *
	 * @param  tree1  First DNA sequence tree.
	 * @param  tree2  Second DNA sequence tree.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>tree1</TT> is null. Thrown if
	 *     <TT>tree2</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if this tree's capacity is less than
	 *     <I>N</I><SUB>1</SUB>&nbsp;+&nbsp;<I>N</I><SUB>2</SUB>&nbsp;+&nbsp;1.
	 */
	public void join
		(DnaSequenceTree tree1,
		 DnaSequenceTree tree2)
		{
		// Verify preconditions.
		int N1 = tree1.myLength;
		int N2 = tree2.myLength;
		if (this.myCapacity < N1+N2+1)
			{
			throw new IllegalArgumentException
				("DnaSequenceTree.join(): Capacity (= "+this.myCapacity+
				 ") too small");
			}
		int C = this.myCapacity;

		// Set up this tree's root node at index 0.
		Node root = this.myNode[0];
		root.parent = -1;
		root.child1 = tree1.myRoot + 1;
		root.child2 = tree2.myRoot + N1 + 1;
		root.seq = null;
		root.brlen = null;

		// Copy <tree1>'s nodes to this tree. Must offset <tree1>'s node
		// indexes by 1.
		for (int i = 0; i < N1; ++ i)
			{
			Node thisnode = this.myNode[i+1];
			Node treenode = tree1.myNode[i];
			thisnode.parent =
				treenode.parent == -1 ? 0 : treenode.parent + 1;
			thisnode.child1 =
				treenode.child1 == -1 ? -1 : treenode.child1 + 1;
			thisnode.child2 =
				treenode.child2 == -1 ? -1 : treenode.child2 + 1;
			thisnode.seq = treenode.seq;
			thisnode.brlen = treenode.brlen;
			}

		// Copy <tree2>'s nodes to this tree. Must offset <tree2>'s node
		// indexes by N1 + 1.
		for (int i = 0; i < N2; ++ i)
			{
			Node thisnode = this.myNode[i+N1+1];
			Node treenode = tree2.myNode[i];
			thisnode.parent =
				treenode.parent == -1 ? 0 : treenode.parent + N1 + 1;
			thisnode.child1 =
				treenode.child1 == -1 ? -1 : treenode.child1 + N1 + 1;
			thisnode.child2 =
				treenode.child2 == -1 ? -1 : treenode.child2 + N1 + 1;
			thisnode.seq = treenode.seq;
			thisnode.brlen = treenode.brlen;
			}

		// Clear any unused nodes in this tree.
		for (int i = N1+N2+1; i < C; ++ i)
			{
			this.myNode[i].clear();
			}

		// Update this tree's length and root.
		this.myLength = N1+N2+1;
		this.myRoot = 0;
		}

	/**
	 * Add the given DNA sequence to this DNA sequence tree. This tree's
	 * capacity is unchanged and must be greater than or equal to this tree's
	 * length + 2. A new node is added to this tree. No DNA sequence is
	 * associated with the new node. The new node's parent is the parent of the
	 * node at index <TT>i</TT>; if the node at index <TT>i</TT> was the root
	 * node, the new node becomes the root node. The new node's first child is
	 * the node at index <TT>i</TT>. The new node's second child is a new tip
	 * node associated with the given DNA sequence.
	 * <P>
	 * Alternatively, if this tree is empty, the <TT>add()</TT> method sets this
	 * tree to have one root node associated with the given DNA sequence. In
	 * this case <TT>i</TT> is ignored.
	 * <P>
	 * <I>Note:</I> This tree contains a reference to (not a copy of)
	 * <TT>seq</TT>.
	 * <P>
	 * <I>Note:</I> This method may alter the index of this tree's root node.
	 *
	 * @param  i    Node index, 0 &le; <TT>i</TT> &le; <I>N</I>&minus;1.
	 * @param  seq  DNA sequence associated with new tip node, or null if no DNA
	 *              sequence is associated.
	 *
	 * @return  Index of new tip node.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if this tree's capacity is less than
	 *     this tree's length + 2.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public int add
		(int i,
		 DnaSequence seq)
		{
		// Case 1: This tree is empty.
		if (myLength == 0)
			{
			// Verify preconditions.
			if (myCapacity < 1)
				{
				throw new IllegalArgumentException
					("DnaSequenceTree.add(): Capacity (= "+myCapacity+
					 ") < 1, illegal");
				}

			// Set up root node with <seq>.
			Node root = this.myNode[0];
			root.parent = -1;
			root.child1 = -1;
			root.child2 = -1;
			root.seq = seq;
			root.brlen = null;

			// Update this tree's length and root.
			myLength = 1;
			myRoot = 0;

			return 0;
			}

		// Case 2: This tree is not empty.
		else
			{
			// Verify preconditions.
			if (myCapacity < myLength+2)
				{
				throw new IllegalArgumentException
					("DnaSequenceTree.add(): Capacity (= "+myCapacity+
					 ") too small");
				}
			if (0 > i || i >= myLength)
				{
				throw new IndexOutOfBoundsException
					("DnaSequenceTree.add(): Index (= "+i+") out of bounds");
				}

			// Set up new interior node.
			int parent_i = myNode[i].parent;
			Node newnode = myNode[myLength];
			newnode.parent = parent_i;
			newnode.child1 = i;
			newnode.child2 = myLength+1;
			newnode.seq = null;
			newnode.brlen = null;

			// Splice new interior node between node <i> and node <i>'s parent.
			// Update this tree's root if necessary.
			if (parent_i == -1)
				{
				myRoot = myLength;
				}
			else if (myNode[parent_i].child1 == i)
				{
				myNode[parent_i].child1 = myLength;
				}
			else
				{
				myNode[parent_i].child2 = myLength;
				}
			myNode[i].parent = myLength;

			// Set up new tip node with <seq>.
			Node newtip = myNode[myLength+1];
			newtip.parent = myLength;
			newtip.child1 = -1;
			newtip.child2 = -1;
			newtip.seq = seq;
			newtip.brlen = null;

			// Update this tree's length.
			myLength += 2;

			return myLength - 1;
			}
		}

	/**
	 * Create a {@linkplain DnaSequenceList} consisting of the DNA sequences
	 * associated with the tip nodes in this tree. The DNA sequences appear in
	 * descending order of the tip nodes' branch lengths. A tip node with no
	 * associated branch length uses a default branch length of 0.
	 * <P>
	 * <I>Note:</I> The returned list contains references to (not copies of)
	 * the DNA sequences in this tree.
	 *
	 * @return  DNA sequence list.
	 */
	public DnaSequenceList toList()
		{
		// Make a list of just the tip nodes.
		ArrayList<Node> nodelist = new ArrayList<Node>();
		int N = myLength;
		for (int i = 0; i < N; ++ i)
			{
			Node node = myNode[i];
			if (node.child1 == -1) nodelist.add (node);
			}

		// Sort the node list.
		Collections.sort (nodelist);
		int M = nodelist.size();

		// Set up DNA sequence list.
		DnaSequenceList seqlist = new DnaSequenceList();
		seqlist.mySequence = new DnaSequence [M];
		int i = 0;
		for (Node node : nodelist)
			{
			seqlist.mySequence[i++] = node.seq;
			}

		return seqlist;
		}

	/**
	 * Returns a string version of this DNA sequence tree. The returned string
	 * is in Newick Standard format, including the branch lengths if any and the
	 * tip nodes' DNA sequence names. For further information about Newick
	 * Standard format, see:
	 * <UL>
	 * <LI>
	 * <A HREF="http://evolution.gs.washington.edu/phylip/newicktree.html" TARGET="_top">http://evolution.gs.washington.edu/phylip/newicktree.html</A>
	 * </UL>
	 *
	 * @return  String version of this tree.
	 */
	public String toString()
		{
		StringBuilder buf = new StringBuilder();
		if (myLength == 0)
			{
			buf.append ('(');
			buf.append (')');
			}
		else if (myLength == 1)
			{
			buf.append ('(');
			toString (buf, 0);
			buf.append (')');
			}
		else
			{
			toString (buf, myRoot);
			}
		buf.append (';');
		return buf.toString();
		}

	private void toString
		(StringBuilder buf,
		 int index)
		{
		Node node = myNode[index];
		if (node.child1 != -1)
			{
			buf.append ('(');
			toString (buf, node.child1);
			buf.append (',');
			toString (buf, node.child2);
			buf.append (')');
			}
		buf.append (nodeName (node));
		if (node.brlen != null)
			{
			buf.append (':');
			buf.append (node.brlen);
			}
		}

	private String nodeName
		(Node node)
		{
		return
			node.seq == null || node.seq.myName == null ?
				"" :
				node.seq.myName.replaceAll ("\\s+", "_");
		}

	}
