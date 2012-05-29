//******************************************************************************
//
// File:    DynamicIntegerSchedule.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.DynamicIntegerSchedule
//
// This Java source file is copyright (C) 2009 by Alan Kaminsky. All rights
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

package edu.rit.pj;

import java.util.concurrent.atomic.AtomicInteger;

import edu.rit.util.Range;

/**
 * Class DynamicIntegerSchedule provides a dynamic schedule object. The loop
 * index is type <TT>int</TT>. The loop iterations are apportioned into chunks
 * of a given size (a given number of iterations per chunk). Each parallel team
 * thread repeatedly performs the next available chunk of iterations until there
 * are no more chunks. The final chunk may be smaller than the given chunk size.
 *
 * @author  Alan Kaminsky
 * @version 17-Nov-2009
 */
class DynamicIntegerSchedule
	extends IntegerSchedule
	{

// Hidden data members.

	// Loop iteration range.
	private Range myLoopRange;

	// Number of iterations already handed out.
	private AtomicInteger N1 = new AtomicInteger();

	// Chunk size.
	private int N2;

// Exported constructors.

	/**
	 * Construct a new dynamic schedule object with a chunk size of 1.
	 */
	public DynamicIntegerSchedule()
		{
		this (1);
		}

	/**
	 * Construct a new dynamic schedule object with the given chunk size.
	 *
	 * @param  theChunkSize  Chunk size.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theChunkSize</TT> is less than 1.
	 */
	public DynamicIntegerSchedule
		(int theChunkSize)
		{
		super();
		if (theChunkSize < 1)
			{
			throw new IllegalArgumentException
				("DynamicIntegerSchedule(): Chunk size = " + theChunkSize +
				 " illegal");
			}
		N2 = theChunkSize;
		}

	/**
	 * Construct a new dynamic schedule object. This constructor is for use by
	 * the <TT>IntegerSchedule.parse()</TT> method. <TT>args</TT> must be an
	 * array of one string, namely the chunk size, an integer &gt;= 1.
	 *
	 * @param  args  Array of argument strings.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>args</TT> is not an array of one
	 *     string. Thrown if the chunk size is less than 1.
	 */
	public DynamicIntegerSchedule
		(String[] args)
		{
		this (getChunkSize (args));
		}

	private static int getChunkSize
		(String[] args)
		{
		if (args.length != 1)
			{
			throw new IllegalArgumentException
				("DynamicIntegerSchedule(): Usage: -Dpj.schedule=dynamic or -Dpj.schedule=\"dynamic(<n>)\"");
			}
		int theChunkSize;
		try
			{
			theChunkSize = Integer.parseInt (args[0]);
			}
		catch (NumberFormatException exc)
			{
			throw new IllegalArgumentException
				("DynamicIntegerSchedule(): Chunk size = " + args[0] +
				 " illegal");
			}
		return theChunkSize;
		}

	/**
	 * Determine if this schedule is a fixed schedule. For a parallel team with
	 * <I>K</I> threads, a fixed schedule partitions the loop index range into
	 * exactly <I>K</I> chunks, one chunk for each thread, each chunk with
	 * predetermined upper and lower bounds.
	 *
	 * @return  True if this is a fixed schedule, false otherwise.
	 */
	public boolean isFixedSchedule()
		{
		return false;
		}

// Hidden operations.

	/**
	 * Start generating chunks of iterations for a parallel for loop using this
	 * schedule.
	 * <P>
	 * The <TT>start()</TT> method is only called by a single thread in the
	 * Parallel Java middleware.
	 *
	 * @param  K             Number of threads in the parallel team.
	 * @param  theLoopRange  Range of iterations for the entire parallel for
	 *                       loop. The stride may be 1 or greater.
	 */
	public void start
		(int K,
		 Range theLoopRange)
		{
		myLoopRange = theLoopRange;
		N1.set (0);
		}

	/**
	 * Obtain the next chunk of iterations for the given thread index. If there
	 * are more iterations, a range object is returned whose lower bound, upper
	 * bound, and stride specify the chunk of iterations to perform. The
	 * returned range object's stride is the same as that given to the
	 * <TT>start()</TT> method. The returned range object's lower bound and
	 * upper bound are contained within the range given to the <TT>start()</TT>
	 * method. If there are no more iterations, null is returned.
	 * <P>
	 * The <TT>next()</TT> method is called by multiple parallel team threads in
	 * the Parallel Java middleware. The <TT>next()</TT> method must be multiple
	 * thread safe.
	 *
	 * @param  theThreadIndex  Thread index in the range 0 .. <I>K</I>-1.
	 *
	 * @return  Chunk of iterations, or null if no more iterations.
	 */
	public Range next
		(int theThreadIndex)
		{
		for (;;)
			{
			int oldN1 = N1.get();
			Range result = myLoopRange.chunk (oldN1, N2);
			int N = result.length();
			if (N == 0) return null;
			int newN1 = oldN1 + N;
			if (N1.compareAndSet (oldN1, newN1)) return result;
			}
		}

	}
