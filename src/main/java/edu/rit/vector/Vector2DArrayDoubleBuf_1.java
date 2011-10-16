//******************************************************************************
//
// File:    Vector2DArrayDoubleBuf_1.java
// Package: edu.rit.vector
// Unit:    Class edu.rit.vector.Vector2DArrayDoubleBuf_1
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

package edu.rit.vector;

import edu.rit.mp.Buf;
import edu.rit.mp.DoubleBuf;

import edu.rit.pj.reduction.Op;

import edu.rit.util.Range;

import java.nio.ByteBuffer;

/**
 * Class Vector2DArrayDoubleBuf_1 provides a buffer for an array of {@linkplain
 * Vector2D} items sent or received using the Message Protocol (MP). The array
 * element stride must be 1. In a message, a vector is treated as two double
 * values (<I>x</I>,<I>y</I>), with each double value represented as eight
 * bytes, most significant byte first.
 * <P>
 * A buffer may be used to send one or more messages at the same time in
 * multiple threads. If a buffer is being used to send a message or messages,
 * the buffer must not be used to receive a message at the same time.
 * <P>
 * A buffer may be used to receive one message at a time. If a buffer is being
 * used to receive a message, the buffer must not be used to receive another
 * message in a different thread, and the buffer must not be used to send a
 * message or messages.
 * <P>
 * A buffer is a conduit for retrieving and storing data in some underlying data
 * structure. If the underlying data structure is multiple thread safe, then one
 * thread can be retrieving or storing data via the buffer at the same time as
 * other threads are accessing the data structure. If the underlying data
 * structure is not multiple thread safe, then other threads must not access the
 * data structure while one thread is retrieving or storing data via the buffer.
 *
 * @author  Alan Kaminsky
 * @version 27-Jun-2007
 */
class Vector2DArrayDoubleBuf_1
	extends Vector2DArrayDoubleBuf
	{

// Hidden constructors.

	/**
	 * Construct a new vector array buffer.
	 *
	 * @param  theArray  Array.
	 * @param  theRange  Range of elements to include. The stride is assumed to
	 *                   be 1.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theArray</TT>'s allocation does
	 *     not include <TT>theRange</TT>.
	 */
	public Vector2DArrayDoubleBuf_1
		(Vector2D[] theArray,
		 Range theRange)
		{
		super (theArray, theRange);
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
	 */
	public double get
		(int i)
		{
		int index = myArrayOffset + (i >>> 1);
		if ((i & 1) == 0)
			{
			return myArray[index].x;
			}
		else
			{
			return myArray[index].y;
			}
		}

	/**
	 * Store the given item in this buffer.
	 * <P>
	 * The <TT>put()</TT> method must not block the calling thread; if it does,
	 * all message I/O in MP will be blocked.
	 *
	 * @param  i     Item index in the range 0 .. <TT>length()</TT>-1.
	 * @param  item  Item to be stored at index <TT>i</TT>.
	 */
	public void put
		(int i,
		 double item)
		{
		int index = myArrayOffset + (i >>> 1);
		if ((i & 1) == 0)
			{
			myArray[index].x = item;
			}
		else
			{
			myArray[index].y = item;
			}
		}

	/**
	 * Create a buffer for performing parallel reduction using the given binary
	 * operation. The results of the reduction are placed into this buffer.
	 * <P>
	 * <I>Note:</I> Class Vector2DArrayBuf_1 does not support reduction. The
	 * <TT>getReductionBuf()</TT> method throws an
	 * UnsupportedOperationException.
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
		int index = myArrayOffset + (i >>> 1);
		int off = i;
		while (off < myLength && buffer.remaining() >= 16)
			{
			Vector2D item = myArray[index];
			buffer.putDouble (item.x);
			buffer.putDouble (item.y);
			++ index;
			off += 2;
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
		int index = myArrayOffset + (i >>> 1);
		int off = i;
		int max = Math.min (i + num, myLength);
		while (off < max && buffer.remaining() >= 16)
			{
			Vector2D item = myArray[index];
			item.x = buffer.getDouble();
			item.y = buffer.getDouble();
			++ index;
			off += 2;
			}
		return off - i;
		}

	}
