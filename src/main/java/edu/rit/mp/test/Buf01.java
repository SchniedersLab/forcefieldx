//******************************************************************************
//
// File:    Buf01.java
// Package: edu.rit.mp.test
// Unit:    Class edu.rit.mp.test.Buf01
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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

package edu.rit.mp.test;

import edu.rit.pj.reduction.Op;

import edu.rit.mp.Buf;
import edu.rit.mp.IntegerBuf;

import java.nio.ByteBuffer;

/**
 * Class Buf01 provides an integer buffer for the {@linkplain Send01} and
 * {@linkplain Receive01} programs.
 *
 * @author  Alan Kaminsky
 * @version 24-Jun-2007
 */
public class Buf01
	extends IntegerBuf
	{

// Exported constructors.

	/**
	 * Construct a new integer buffer.
	 *
	 * @param  theLength     Number of items.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theLength</TT> &lt; 0.
	 */
	public Buf01
		(int theLength)
		{
		super (theLength);
		}

// Exported operations.

	/**
	 * Obtain the given item from this buffer.
	 * <P>
	 * The <TT>get()</TT> method must not block the calling thread; if it does,
	 * all message I/O in MP will be blocked.
	 *
	 * @param  i  Item index in the range 0 .. <TT>length()</TT>-1.
	 *
	 * @return  Item at index <TT>i</TT>.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public int get
		(int i)
		{
		return i;
		}

	/**
	 * Store the given item in this buffer.
	 * <P>
	 * The <TT>put()</TT> method must not block the calling thread; if it does,
	 * all message I/O in MP will be blocked.
	 *
	 * @param  i     Item index in the range 0 .. <TT>length()</TT>-1.
	 * @param  item  Item to be stored at index <TT>i</TT>.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public void put
		(int i,
		 int item)
		{
		if (i != item)
			{
			System.out.print ("Item ");
			System.out.print (i);
			System.out.print (", value = ");
			System.out.print (item);
			System.out.println();
			}
		}

	/**
	 * Create a buffer for performing parallel reduction using the given binary
	 * operation. The results of the reduction are placed into this buffer.
	 * <P>
	 * Operations performed on the returned reduction buffer have the same
	 * effect as operations performed on this buffer, except whenever a source
	 * item <I>S</I> is put into a destination item <I>D</I> in this buffer,
	 * <I>D</I> is set to <I>D op S</I>, that is, the reduction of <I>D</I> and
	 * <I>S</I> using the given binary operation (rather than just setting
	 * <I>D</I> to <I>S</I>).
	 *
	 * @param  op  Binary operation.
	 *
	 * @exception  ClassCastException
	 *     (unchecked exception) Thrown if this buffer's element data type and
	 *     the given binary operation's argument data type are not the same.
	 */
	public Buf getReductionBuf
		(Op op)
		{
		throw new UnsupportedOperationException();
		}

// Hidden operations.

	/**
	 * Send as many items as possible from this buffer to the given byte
	 * buffer.
	 * <P>
	 * The <TT>sendItems()</TT> method must not block the calling thread; if it
	 * does, all message I/O in MP will be blocked.
	 *
	 * @param  i       Index of first item to send, in the range 0 ..
	 *                 <TT>length</TT>-1.
	 * @param  buffer  Byte buffer.
	 *
	 * @return  Number of items sent.
	 */
	protected int sendItems
		(int i,
		 ByteBuffer buffer)
		{
		int off = i;
		int len = length();
		while (off < len && buffer.remaining() >= 4)
			{
			buffer.putInt (off);
			++ off;
			}
		return off - i;
		}

	/**
	 * Receive as many items as possible from the given byte buffer to this
	 * buffer.
	 * <P>
	 * The <TT>receiveItems()</TT> method must not block the calling thread; if
	 * it does, all message I/O in MP will be blocked.
	 *
	 * @param  i       Index of first item to receive, in the range 0 ..
	 *                 <TT>length</TT>-1.
	 * @param  num     Maximum number of items to receive.
	 * @param  buffer  Byte buffer.
	 *
	 * @return  Number of items received.
	 */
	protected int receiveItems
		(int i,
		 int num,
		 ByteBuffer buffer)
		{
		int off = i;
		int max = Math.min (i + num, length());
		while (off < max && buffer.remaining() >= 4)
			{
			put (off, buffer.getInt());
			++ off;
			}
		return off - i;
		}

	}
