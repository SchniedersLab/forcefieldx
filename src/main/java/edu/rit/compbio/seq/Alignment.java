//******************************************************************************
//
// File:    Alignment.java
// Package: edu.rit.compbio.seq
// Unit:    Class edu.rit.compbio.seq.Alignment
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

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class Alignment encapsulates the result of a sequence alignment between a
 * query sequence and a subject sequence. Alignment objects are returned by the
 * <TT>align()</TT> method of class {@linkplain ProteinLocalAlignment}.
 * <P>
 * An alignment object does not record the actual query sequence and subject
 * sequence that were aligned. The alignment object does record an ID associated
 * with the query sequence and an ID associated with the subject sequence. For
 * example, the ID could be the index of a sequence in a {@linkplain
 * ProteinDatabase}.
 * <P>
 * Class Alignment implements interface Comparable and so has a natural
 * ordering. The comparison order depends on the alignment scores and the
 * subject sequence IDs. An alignment with a higher score comes before an
 * alignment with a lower score. If the scores are equal, an alignment with a
 * lower subject sequence ID comes before an alignment with a higher subject
 * sequence ID.
 *
 * @author  Alan Kaminsky
 * @version 03-Jul-2008
 */
public class Alignment
	implements Externalizable, Comparable<Alignment>
	{

// Exported constants.

	/**
	 * State of an alignment position: character in query sequence aligned with
	 * character in subject sequence.
	 */
	public static final int QUERY_ALIGNED_WITH_SUBJECT = 0;

	/**
	 * State of an alignment position: character in query sequence aligned with
	 * gap in subject sequence.
	 */
	public static final int QUERY_ALIGNED_WITH_GAP = 1;

	/**
	 * State of an alignment position: character in subject sequence aligned
	 * with gap in query sequence.
	 */
	public static final int SUBJECT_ALIGNED_WITH_GAP = 2;

// Hidden data members.

	private static final long serialVersionUID = 8656834485638780975L;

	// Sequence IDs.
	long myQueryId;
	long mySubjectId;

	// Query and subject sequence length.
	int myQueryLength;
	int mySubjectLength;

	// Alignment score.
	int myScore;

	// Location of alignment start point.
	int myQueryStart;
	int mySubjectStart;

	// Location of alignment finish point.
	int myQueryFinish;
	int mySubjectFinish;

	// Traceback information (in reverse order).
	byte[] myTraceback;

// Exported constructors.

	/**
	 * Construct a new, uninitialized alignment object. This constructor is for
	 * use only by object deserialization.
	 */
	public Alignment()
		{
		}

// Exported operations.

	/**
	 * Get the ID associated with the query sequence. This is an arbitrary
	 * integer; for example, it could be the index of a sequence in a
	 * {@linkplain ProteinDatabase}.
	 *
	 * @return  Query sequence ID.
	 */
	public long getQueryId()
		{
		return myQueryId;
		}

	/**
	 * Get the ID associated with the subject sequence. This is an arbitrary
	 * integer; for example, it could be the index of a sequence in a
	 * {@linkplain ProteinDatabase}.
	 *
	 * @return  Subject sequence ID.
	 */
	public long getSubjectId()
		{
		return mySubjectId;
		}

	/**
	 * Get the length of the query sequence.
	 *
	 * @return  Query sequence length.
	 */
	public int getQueryLength()
		{
		return myQueryLength;
		}

	/**
	 * Get the length of the subject sequence.
	 *
	 * @return  Subject sequence length.
	 */
	public int getSubjectLength()
		{
		return mySubjectLength;
		}

	/**
	 * Get the index of the first aligned character in the query sequence. The
	 * index is in the range 1 .. <I>M</I>, where <I>M</I> is the length of the
	 * query sequence. If a local alignment was not found, 0 is returned.
	 *
	 * @return  Query sequence start index.
	 */
	public int getQueryStart()
		{
		return myQueryStart;
		}

	/**
	 * Get the index of the last aligned character in the query sequence. The
	 * index is in the range 1 .. <I>M</I>, where <I>M</I> is the length of the
	 * query sequence. If a local alignment was not found, 0 is returned.
	 *
	 * @return  Query sequence finish index.
	 */
	public int getQueryFinish()
		{
		return myQueryFinish;
		}

	/**
	 * Get the index of the first aligned character in the subject sequence.
	 * The index is in the range 1 .. <I>N</I>, where <I>N</I> is the length of
	 * the subject sequence. If a local alignment was not found, 0 is returned.
	 *
	 * @return  Subject sequence start index.
	 */
	public int getSubjectStart()
		{
		return mySubjectStart;
		}

	/**
	 * Get the index of the last aligned character in the subject sequence. The
	 * index is in the range 1 .. <I>N</I>, where <I>N</I> is the length of the
	 * subject sequence. If a local alignment was not found, 0 is returned.
	 *
	 * @return  Subject sequence finish index.
	 */
	public int getSubjectFinish()
		{
		return mySubjectFinish;
		}

	/**
	 * Get the number of positions in the alignment. If a local alignment was
	 * not found, 0 is returned.
	 *
	 * @return  Alignment length.
	 */
	public int getAlignmentLength()
		{
		return myTraceback.length;
		}

	/**
	 * Get the state of the given position in the alignment. The index
	 * <TT>i</TT> must be in the range 0 .. <I>L</I>&minus;1, where <I>L</I> is
	 * the alignment length. The state of position <TT>i</TT> in the alignment
	 * is returned, one of the following:
	 * <UL>
	 * <LI>{@link #QUERY_ALIGNED_WITH_SUBJECT}
	 * <LI>{@link #QUERY_ALIGNED_WITH_GAP}
	 * <LI>{@link #SUBJECT_ALIGNED_WITH_GAP}
	 * </UL>
	 *
	 * @param  i  Position in the alignment.
	 *
	 * @return  State of position <TT>i</TT> in the alignment.
	 */
	public int getAlignment
		(int i)
		{
		return myTraceback[myTraceback.length-1-i];
		}

	/**
	 * Compare this alignment object to the given alignment object. The
	 * comparison order depends on the alignment scores and the subject sequence
	 * IDs. An alignment with a higher score comes before an alignment with a
	 * lower score. If the scores are equal, an alignment with a lower subject
	 * sequence ID comes before an alignment with a higher subject sequence ID.
	 *
	 * @param  alignment  Alignment to compare to.
	 *
	 * @return  An integer less than, equal to, or greater than 0 if this
	 *          alignment comes before, is the same as, or comes after the given
	 *          alignment, respectively.
	 */
	public int compareTo
		(Alignment alignment)
		{
		if (this.myScore > alignment.myScore) return -1;
		else if (this.myScore < alignment.myScore) return 1;
		else if (this.mySubjectId < alignment.mySubjectId) return -1;
		else if (this.mySubjectId > alignment.mySubjectId) return 1;
		else return 0;
		}

	/**
	 * Returns a string version of this alignment object.
	 *
	 * @return  String version.
	 */
	public String toString()
		{
		StringBuilder b = new StringBuilder();
		b.append ("Alignment(qid=");
		b.append (myQueryId);
		b.append (",sid=");
		b.append (mySubjectId);
		b.append (",qlen=");
		b.append (myQueryLength);
		b.append (",slen=");
		b.append (mySubjectLength);
		b.append (",score=");
		b.append (myScore);
		b.append (",qstart=");
		b.append (myQueryStart);
		b.append (",sstart=");
		b.append (mySubjectStart);
		b.append (",qfin=");
		b.append (myQueryFinish);
		b.append (",sfin=");
		b.append (mySubjectFinish);
		b.append (")");
		return b.toString();
		}

	/**
	 * Write this alignment object to the given object output stream.
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
		out.writeLong (myQueryId);
		out.writeLong (mySubjectId);
		out.writeInt (myQueryLength);
		out.writeInt (mySubjectLength);
		out.writeInt (myScore);
		out.writeInt (myQueryStart);
		out.writeInt (mySubjectStart);
		out.writeInt (myQueryFinish);
		out.writeInt (mySubjectFinish);
		int n = myTraceback.length;
		out.writeInt (n);
		for (int i = 0; i < n; ++ i)
			{
			out.writeByte (myTraceback[i]);
			}
		}

	/**
	 * Read this alignment object from the given object input stream.
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
		myQueryId = in.readLong();
		mySubjectId = in.readLong();
		myQueryLength = in.readInt();
		mySubjectLength = in.readInt();
		myScore = in.readInt();
		myQueryStart = in.readInt();
		mySubjectStart = in.readInt();
		myQueryFinish = in.readInt();
		mySubjectFinish = in.readInt();
		int n = in.readInt();
		myTraceback = new byte [n];
		for (int i = 0; i < n; ++ i)
			{
			myTraceback[i] = in.readByte();
			}
		}

	}
