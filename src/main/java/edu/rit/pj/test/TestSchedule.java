//******************************************************************************
//
// File:    TestSchedule.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.TestSchedule
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

package edu.rit.pj.test;

import edu.rit.pj.IntegerSchedule;

import edu.rit.util.Range;

import java.util.concurrent.atomic.AtomicLong;

/**
 * Class TestSchedule provides a schedule object for testing class {@linkplain
 * edu.rit.pj.IntegerSchedule}.
 *
 * @author  Alan Kaminsky
 * @version 17-Nov-2009
 */
public class TestSchedule
	extends IntegerSchedule
	{

// Hidden data members.

	// Loop iteration range.
	private Range myLoopRange;

	// Number of iterations already handed out, N1 (upper 32 bits) plus chunk
	// size, N2 (lower 32 bits).
	private AtomicLong N1N2 = new AtomicLong();

	// Initial chunk size.
	private int initN2;

	// Chunk size increment.
	private int INCR;

// Exported constructors.

	/**
	 * Construct a new test schedule object with an initial chunk size of 1 and
	 * a chunk size increment of 0.
	 */
	public TestSchedule()
		{
		this (1, 0);
		}

	/**
	 * Construct a new test schedule object with the given initial chunk size
	 * and chunk size increment.
	 *
	 * @param  theInitialChunkSize    Initial chunk size.
	 * @param  theChunkSizeIncrement  Chunk size increment.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theInitialChunkSize</TT> is less
	 *     than 1. Thrown if <TT>theChunkSizeIncrement</TT> is less than 0.
	 */
	public TestSchedule
		(int theInitialChunkSize,
		 int theChunkSizeIncrement)
		{
		super();
		if (theInitialChunkSize < 1)
			{
			throw new IllegalArgumentException
				("TestSchedule(): Initial chunk size = " +
				 theInitialChunkSize + " illegal");
			}
		if (theChunkSizeIncrement < 0)
			{
			throw new IllegalArgumentException
				("TestSchedule(): Chunk size increment = " +
				 theChunkSizeIncrement + " illegal");
			}
		initN2 = theInitialChunkSize;
		INCR = theChunkSizeIncrement;
		}

	/**
	 * Construct a new test schedule object. This constructor is for use by the
	 * <TT>IntegerSchedule.parse()</TT> method. <TT>args</TT> must be an array
	 * of two strings, namely the initial chunk size, an integer &gt;= 1, and
	 * the chunk size increment, an integer &gt;= 0.
	 *
	 * @param  args  Array of argument strings.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>args</TT> is not an array of two
	 *     strings. Thrown if the initial chunk size is less than 1. Thrown if
	 *     the chunk size increment is less than 0.
	 */
	public TestSchedule
		(String[] args)
		{
		this (getInitialChunkSize (args), getChunkSizeIncrement (args));
		}

	private static int getInitialChunkSize
		(String[] args)
		{
		if (args.length != 2)
			{
			throw new IllegalArgumentException
				("TestSchedule(): Usage: -Dpj.schedule=edu.rit.pj.test.TestSchedule or -Dpj.schedule=\"edu.rit.pj.test.TestSchedule(<n>,<incr>)\"");
			}
		int theInitialChunkSize;
		try
			{
			theInitialChunkSize = Integer.parseInt (args[0]);
			}
		catch (NumberFormatException exc)
			{
			throw new IllegalArgumentException
				("TestSchedule(): Initial chunk size = " + args[0] +
				 " illegal");
			}
		return theInitialChunkSize;
		}

	private static int getChunkSizeIncrement
		(String[] args)
		{
		if (args.length != 2)
			{
			throw new IllegalArgumentException
				("TestSchedule(): Usage: -Dpj.schedule=edu.rit.pj.test.TestSchedule or -Dpj.schedule=\"edu.rit.pj.test.TestSchedule(<n>,<incr>)\"");
			}
		int theChunkSizeIncrement;
		try
			{
			theChunkSizeIncrement = Integer.parseInt (args[1]);
			}
		catch (NumberFormatException exc)
			{
			throw new IllegalArgumentException
				("TestSchedule(): Chunk size increment = " + args[0] +
				 " illegal");
			}
		return theChunkSizeIncrement;
		}

// Exported operations.

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
	 * Start a parallel for loop using this schedule.
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
		N1N2.set (initN2 & 0x00000000FFFFFFFFL);
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
			long oldN1N2 = N1N2.get();
			int oldN1 = (int) (oldN1N2 >>> 32);
			int oldN2 = (int) oldN1N2;
			Range result = myLoopRange.chunk (oldN1, oldN2);
			int N = result.length();
			if (N == 0) return null;
			int newN1 = oldN1 + N;
			int newN2 = oldN2 + INCR;
			long newN1N2 =
				((newN1 & 0x00000000FFFFFFFFL) << 32) |
				 (newN2 & 0x00000000FFFFFFFFL);
			if (N1N2.compareAndSet (oldN1N2, newN1N2)) return result;
			}
		}

	}
