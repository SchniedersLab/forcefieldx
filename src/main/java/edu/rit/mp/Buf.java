//******************************************************************************
//
// File:    Buf.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.Buf
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

package edu.rit.mp;

import edu.rit.pj.reduction.Op;

import java.io.IOException;

import java.nio.ByteBuffer;

/**
 * Class Buf is the abstract base class for a buffer of items sent or received
 * using the Message Protocol (MP).
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
 * @version 03-May-2008
 */
public abstract class Buf
	{

// Hidden data members.

	/**
	 * Number of items in this buffer.
	 */
	protected final int myLength;

	// Type and length sent or received in a message.
	byte myMessageType;
	int myMessageLength;

// Hidden constructors.

	/**
	 * Construct a new buffer.
	 *
	 * @param  theType       Item type.
	 * @param  theLength     Number of items.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theLength</TT> &lt; 0.
	 */
	Buf
		(byte theType,
		 int theLength)
		{
		if (theLength < 0)
			{
			throw new IllegalArgumentException
				("edu.rit.mp.Buf(): theLength = " + theLength + " illegal");
			}
		myLength = theLength;
		myMessageType = theType;
		myMessageLength = theLength;
		}

// Exported operations.

	/**
	 * Obtain the number of items in this buffer.
	 *
	 * @return  Number of items.
	 */
	public final int length()
		{
		return myLength;
		}

	/**
	 * Copy items from the given buffer to this buffer. The number of items
	 * copied is this buffer's length or <TT>theSrc</TT>'s length, whichever is
	 * smaller. If <TT>theSrc</TT> is this buffer, the <TT>copy()</TT> method
	 * does nothing.
	 *
	 * @param  theSrc  Source of items to copy into this buffer.
	 *
	 * @exception  ClassCastException
	 *     (unchecked exception) Thrown if <TT>theSrc</TT>'s item data type is
	 *     not the same as this buffer's item data type.
	 */
	public abstract void copy
		(Buf theSrc);

	/**
	 * Fill this buffer with the given item. The <TT>item</TT> is assigned to
	 * each element in this buffer.
	 * <P>
	 * If this buffer's item data type is a primitive type, the <TT>item</TT>
	 * must be an instance of the corresponding primitive wrapper class -- class
	 * Integer for type <TT>int</TT>, class Double for type <TT>double</TT>, and
	 * so on. If the <TT>item</TT> is null, the item data type's default initial
	 * value is assigned to each element in this buffer.
	 * <P>
	 * If this buffer's item data type is a nonprimitive type, the <TT>item</TT>
	 * must be an instance of the item class or a subclass thereof. The
	 * <TT>item</TT> may be null. Note that since <TT>item</TT> is
	 * <I>assigned</I> to every buffer element, every buffer element ends up
	 * referring to the same <TT>item</TT>.
	 *
	 * @param  item  Item.
	 *
	 * @exception  ClassCastException
	 *     (unchecked exception) Thrown if the <TT>item</TT>'s data type is not
	 *     the same as this buffer's item data type.
	 */
	public abstract void fill
		(Object item);

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
	public abstract Buf getReductionBuf
		(Op op);

	/**
	 * Create a temporary buffer with the same type of items and the same length
	 * as this buffer. The new buffer items are stored in a newly created array,
	 * separate from the storage for this buffer's items.
	 */
	public abstract Buf getTemporaryBuf();

// Hidden operations.

	/**
	 * Called by the I/O thread before sending message items using this buffer.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	void preSend()
		throws IOException
		{
		}

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
	protected abstract int sendItems
		(int i,
		 ByteBuffer buffer);

	/**
	 * Called by the I/O thread after sending message items using this buffer.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	void postSend()
		throws IOException
		{
		}

	/**
	 * Called by the I/O thread before receiving message items using this
	 * buffer.
	 *
	 * @param  theReadLength  Actual number of items in message.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	void preReceive
		(int theReadLength)
		throws IOException
		{
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
	protected abstract int receiveItems
		(int i,
		 int num,
		 ByteBuffer buffer);

	/**
	 * Skip as many items as possible from the given byte buffer.
	 *
	 * @param  num     Number of items to skip.
	 * @param  buffer  Buffer.
	 *
	 * @return  Number of items actually skipped.
	 */
	abstract int skipItems
		(int num,
		 ByteBuffer buffer);

	/**
	 * Called by the I/O thread after receiving message items using this
	 * buffer.
	 *
	 * @param  theStatus       Status object that will be returned for the
	 *                         message; its contents may be altered if
	 *                         necessary.
	 * @param  theClassLoader  Alternate class loader to be used when receiving
	 *                         objects, or null.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	void postReceive
		(Status theStatus,
		 ClassLoader theClassLoader)
		throws IOException
		{
		}

	}
