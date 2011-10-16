//******************************************************************************
//
// File:    DnaSequence.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.DnaSequence
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
 * Class DnaSequence encapsulates a DNA sequence. The DNA sequence consists of a
 * sequence of <B>sites</B>. Each site has a <B>state,</B> which is a set of
 * <B>bases</B>. The four bases are adenine, cytosine, guanine, and thymine. For
 * textual I/O, each state is represented by a single character as follows:
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR><TD><I>Char.</I></TD><TD WIDTH=20> </TD>
 * <TD><I>Meaning</I></TD><TD WIDTH=20> </TD>
 * <TD><I>Set</I></TD></TR>
 * <TR><TD>A</TD><TD WIDTH=20> </TD>
 * <TD>Adenine</TD><TD WIDTH=20> </TD>
 * <TD>(A)</TD></TR>
 * <TR><TD>C</TD><TD WIDTH=20> </TD>
 * <TD>Cytosine</TD><TD WIDTH=20> </TD>
 * <TD>(C)</TD></TR>
 * <TR><TD>G</TD><TD WIDTH=20> </TD>
 * <TD>Guanine</TD><TD WIDTH=20> </TD>
 * <TD>(G)</TD></TR>
 * <TR><TD>T</TD><TD WIDTH=20> </TD>
 * <TD>Thymine</TD><TD WIDTH=20> </TD>
 * <TD>(T)</TD></TR>
 * <TR><TD>Y</TD><TD WIDTH=20> </TD>
 * <TD>pYrimidine</TD><TD WIDTH=20> </TD>
 * <TD>(C or T)</TD></TR>
 * <TR><TD>R</TD><TD WIDTH=20> </TD>
 * <TD>puRine</TD><TD WIDTH=20> </TD>
 * <TD>(A or G)</TD></TR>
 * <TR><TD>W</TD><TD WIDTH=20> </TD>
 * <TD>"Weak"</TD><TD WIDTH=20> </TD>
 * <TD>(A or T)</TD></TR>
 * <TR><TD>S</TD><TD WIDTH=20> </TD>
 * <TD>"Strong"</TD><TD WIDTH=20> </TD>
 * <TD>(C or G)</TD></TR>
 * <TR><TD>K</TD><TD WIDTH=20> </TD>
 * <TD>"Keto"</TD><TD WIDTH=20> </TD>
 * <TD>(G or T)</TD></TR>
 * <TR><TD>M</TD><TD WIDTH=20> </TD>
 * <TD>"aMino"</TD><TD WIDTH=20> </TD>
 * <TD>(A or C)</TD></TR>
 * <TR><TD>B</TD><TD WIDTH=20> </TD>
 * <TD>not A</TD><TD WIDTH=20> </TD>
 * <TD>(C or G or T)</TD></TR>
 * <TR><TD>D</TD><TD WIDTH=20> </TD>
 * <TD>not C</TD><TD WIDTH=20> </TD>
 * <TD>(A or G or T)</TD></TR>
 * <TR><TD>H</TD><TD WIDTH=20> </TD>
 * <TD>not G</TD><TD WIDTH=20> </TD>
 * <TD>(A or C or T)</TD></TR>
 * <TR><TD>V</TD><TD WIDTH=20> </TD>
 * <TD>not T</TD><TD WIDTH=20> </TD>
 * <TD>(A or C or G)</TD></TR>
 * <TR><TD>X</TD><TD WIDTH=20> </TD>
 * <TD>unknown</TD><TD WIDTH=20> </TD>
 * <TD>(A or C or G or T)</TD></TR>
 * <TR><TD>-</TD><TD WIDTH=20> </TD>
 * <TD>deletion</TD><TD WIDTH=20> </TD>
 * <TD>()</TD></TR>
 * </TABLE>
 * <P>
 * The DNA sequence has an associated <B>score,</B> an integer. The score can be
 * set to anything and later retrieved.
 * <P>
 * The DNA sequence has a <B>name,</B> a string. The name can be set to anything
 * and later retrieved.
 *
 * @author  Alan Kaminsky
 * @version 23-Jul-2008
 */
public class DnaSequence
	{

// Hidden constants.

	// Amount of extra padding in byte array.
	private static final int PAD = 128;

	// Mapping from the state of a site to the corresponding output character.
	// A=1, C=2, G=4, T=8.
	static final char[] state2char = new char[]
		{/*----*/ '-',
		 /*---A*/ 'A',
		 /*--C-*/ 'C',
		 /*--CA*/ 'M',
		 /*-G--*/ 'G',
		 /*-G-A*/ 'R',
		 /*-GC-*/ 'S',
		 /*-GCA*/ 'V',
		 /*T---*/ 'T',
		 /*T--A*/ 'W',
		 /*T-C-*/ 'Y',
		 /*T-CA*/ 'H',
		 /*TG--*/ 'K',
		 /*TG-A*/ 'D',
		 /*TGC-*/ 'B',
		 /*TGCA*/ 'X'};

	// Mapping from the state of a site to the number of bits turned on at that
	// site. A=1, C=2, G=4, T=8.
	static final int[] state2bitCount = new int[]
		{/*----*/ 0,
		 /*---A*/ 1,
		 /*--C-*/ 1,
		 /*--CA*/ 2,
		 /*-G--*/ 1,
		 /*-G-A*/ 2,
		 /*-GC-*/ 2,
		 /*-GCA*/ 3,
		 /*T---*/ 1,
		 /*T--A*/ 2,
		 /*T-C-*/ 2,
		 /*T-CA*/ 3,
		 /*TG--*/ 2,
		 /*TG-A*/ 3,
		 /*TGC-*/ 3,
		 /*TGCA*/ 4};

// Hidden data members.

	// Sequence data. Each site's set of bases is stored as a bitmap in one
	// byte. A=1, C=2, G=4, T=8. PAD bytes of padding are added to avert cache
	// interference.
	byte[] mySites;

	// Length.
	int myLength;

	// Score.
	int myScore;

	// Name.
	String myName;

	// 128 bytes of extra padding to avert cache interference.
	private transient long p0, p1, p2, p3, p4, p5, p6, p7;
	private transient long p8, p9, pa, pb, pc, pd, pe, pf;

// Exported constructors.

	/**
	 * Construct a new zero-length DNA sequence. The score is initially 0. The
	 * name is initially null.
	 */
	public DnaSequence()
		{
		this (0, 0, null);
		}

	/**
	 * Construct a new DNA sequence with the given length. The score is
	 * initially 0. The name is initially null.
	 *
	 * @param  N  Length (number of sites).
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>N</TT> &lt; 0.
	 */
	public DnaSequence
		(int N)
		{
		this (N, 0, null);
		}

	/**
	 * Construct a new DNA sequence with the given length and score. The name is
	 * initially null.
	 *
	 * @param  N      Length (number of sites).
	 * @param  score  Score.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>N</TT> &lt; 0.
	 */
	public DnaSequence
		(int N,
		 int score)
		{
		this (N, score, null);
		}

	/**
	 * Construct a new DNA sequence with the given length, score, and name.
	 *
	 * @param  N      Length (number of sites).
	 * @param  score  Score.
	 * @param  name   Name. May be null.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>N</TT> &lt; 0.
	 */
	public DnaSequence
		(int N,
		 int score,
		 String name)
		{
		if (N < 0)
			{
			throw new IllegalArgumentException
				("DnaSequence(): N (= "+N+") < 0, illegal)");
			}
		this.mySites = new byte [N+PAD];
		this.myLength = N;
		this.myScore = score;
		this.myName = name;
		}

	/**
	 * Construct a new DNA sequence that is a copy of the given DNA sequence.
	 *
	 * @param  seq  DNA sequence to copy.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>seq</TT> is null.
	 */
	public DnaSequence
		(DnaSequence seq)
		{
		this.mySites = (byte[]) seq.mySites.clone();
		this.myLength = seq.myLength;
		this.myScore = seq.myScore;
		this.myName = seq.myName;
		}

// Exported operations.

	/**
	 * Get this DNA sequence's length.
	 *
	 * @return  Length (number of sites).
	 */
	public int length()
		{
		return myLength;
		}

	/**
	 * Get this DNA sequence's score.
	 *
	 * @return  Score.
	 */
	public int score()
		{
		return myScore;
		}

	/**
	 * Set this DNA sequence's score.
	 *
	 * @param  score  Score.
	 */
	public void score
		(int score)
		{
		myScore = score;
		}

	/**
	 * Get this DNA sequence's name.
	 *
	 * @return  Name. May be null.
	 */
	public String name()
		{
		return myName;
		}

	/**
	 * Set this DNA sequence's name.
	 *
	 * @param  name  Name. May be null.
	 */
	public void name
		(String name)
		{
		myName = name;
		}

	/**
	 * Make this DNA sequence's sites be the same as the given DNA sequence. It
	 * is assumed that this DNA sequence and the given DNA sequence are the same
	 * length. This DNA sequence's score and name are unchanged.
	 *
	 * @param  seq  DNA sequence to copy.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>seq</TT> is null.
	 */
	public void copySites
		(DnaSequence seq)
		{
		System.arraycopy (seq.mySites, 0, this.mySites, 0, myLength);
		}

	/**
	 * Compute the distance between this DNA sequence and the given DNA
	 * sequence. It is assumed that this DNA sequence and the given DNA sequence
	 * are the same length. The distance is the number of differing sites
	 * between the two sequences (the Hamming distance).
	 *
	 * @param  seq  DNA sequence.
	 *
	 * @return  Distance.
	 */
	public double distance
		(DnaSequence seq)
		{
		byte[] site1 = this.mySites;
		byte[] site2 = seq.mySites;
		int diff = 0;
		int N = myLength;
		for (int i = 0; i < N; ++ i)
			{
			if (site1[i] != site2[i]) ++ diff;
			}
		return diff;
		}

	/**
	 * Make this DNA sequence be the ancestor of the two given DNA sequences in
	 * the Fitch parsimony score algorithm. This DNA sequence's sites are set
	 * based on <TT>seq1</TT>'s and <TT>seq2</TT>'s sites. It is assumed that
	 * this DNA sequence and the given DNA sequences are the same length. This
	 * DNA sequence's score is set to the sum of <TT>seq1</TT>'s score,
	 * <TT>seq2</TT>'s score, and the number of state changes at the ancestor.
	 * This DNA sequence's name is unchanged.
	 *
	 * @param  seq1  First child DNA sequence.
	 * @param  seq2  Second child DNA sequence.
	 */
	public void setFitchAncestor
		(DnaSequence seq1,
		 DnaSequence seq2)
		{
		// Get references to sites.
		byte[] ancestor = this.mySites;
		byte[] descendent1 = seq1.mySites;
		byte[] descendent2 = seq2.mySites;
		int N = myLength;

		// Process all sites. Count state changes.
		int nChanges = 0;
		for (int i = 0; i < N; ++ i)
			{
			// Compute intersection of states.
			int state1 = descendent1[i];
			int state2 = descendent2[i];
			int state3 = state1 & state2;

			// If intersection is not empty, record intersection, otherwise
			// record union and note one state change.
			if (state3 == 0)
				{
				state3 = state1 | state2;
				++ nChanges;
				}

			// Update site.
			ancestor[i] = (byte) state3;
			}

		// Record number of state changes.
		this.myScore = seq1.myScore + seq2.myScore + nChanges;
		}

	/**
	 * Returns a string version of this DNA sequence. The string consists of
	 * just the sequence of states (the score and name are not included).
	 */
	public String toString()
		{
		StringBuilder buf = new StringBuilder();
		byte[] site = mySites;
		int N = myLength;
		for (int i = 0; i < N; ++ i)
			{
			buf.append (state2char [site[i]]);
			}
		return buf.toString();
		}

	}
