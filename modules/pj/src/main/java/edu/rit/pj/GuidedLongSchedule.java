//******************************************************************************
//
// File:    GuidedLongSchedule.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.GuidedLongSchedule
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

import edu.rit.util.LongRange;

import java.util.concurrent.atomic.AtomicLong;

/**
 * Class GuidedLongSchedule provides a self-guided schedule object. The loop
 * index is type <TT>long</TT>. The loop iterations are apportioned into chunks
 * of exponentially diminishing sizes. Each successive chunk's size is half the
 * remaining number of iterations divided by the number of threads in the
 * parallel team. However, each chunk is at least a given minimum size (a given
 * minimum number of iterations per chunk). Each parallel team thread repeatedly
 * performs the next available chunk of iterations until there are no more
 * chunks. The final chunk may be smaller than the given minimum chunk size.
 *
 * @author Alan Kaminsky
 * @version 18-Nov-2009
 */
class GuidedLongSchedule
        extends LongSchedule {

// Hidden data members.
    // Twice the number of parallel team threads.
    private int two_K;

    // Loop iteration range.
    private LongRange myLoopRange;

    // Number of iterations.
    private long myLoopRangeLength;

    // Number of iterations already handed out.
    private AtomicLong N1 = new AtomicLong();

    // Minimum chunk size.
    private long N2;

// Exported constructors.
    /**
     * Construct a new self-guided schedule object with a minimum chunk size of
     * 1.
     *
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <TT>theChunkSize</TT> is less than 1.
     */
    public GuidedLongSchedule() {
        this(1);
    }

    /**
     * Construct a new self-guided schedule object.
     *
     * @param theChunkSize Minimum chunk size.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <TT>theChunkSize</TT> is less than 1.
     */
    public GuidedLongSchedule(long theChunkSize) {
        super();
        if (theChunkSize < 1) {
            throw new IllegalArgumentException("GuidedLongSchedule(): Minimum chunk size = "
                    + theChunkSize + " illegal");
        }
        N2 = theChunkSize;
    }

    /**
     * Construct a new self-guided schedule object. This constructor is for use
     * by the <TT>LongSchedule.parse()</TT> method. <TT>args</TT> must be an
     * array of one string, namely the minimum chunk size, an integer &gt;= 1.
     *
     * @param args Array of argument strings.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <TT>args</TT> is not an array of one string. Thrown if the minimum chunk
     * size is less than 1.
     */
    public GuidedLongSchedule(String[] args) {
        this(getChunkSize(args));
    }

    private static long getChunkSize(String[] args) {
        if (args.length != 1) {
            throw new IllegalArgumentException("GuidedLongSchedule(): Usage: -Dpj.schedule=guided or -Dpj.schedule=\"guided(<n>)\"");
        }
        long theChunkSize;
        try {
            theChunkSize = Long.parseLong(args[0]);
        } catch (NumberFormatException exc) {
            throw new IllegalArgumentException("GuidedLongSchedule(): Chunk size = " + args[0]
                    + " illegal");
        }
        return theChunkSize;
    }

// Exported operations.
    /**
     * Determine if this schedule is a fixed schedule. For a parallel team with
     * <I>K</I> threads, a fixed schedule partitions the loop index range into
     * exactly <I>K</I> chunks, one chunk for each thread, each chunk with
     * predetermined upper and lower bounds.
     *
     * @return True if this is a fixed schedule, false otherwise.
     */
    public boolean isFixedSchedule() {
        return false;
    }

// Hidden operations.
    /**
     * {@inheritDoc}
     *
     * Start generating chunks of iterations for a parallel for loop using this
     * schedule.
     * <P>
     * The <TT>start()</TT> method is only called by a single thread in the
     * Parallel Java middleware.
     */
    public void start(int K,
            LongRange theLoopRange) {
        two_K = 2 * K;
        myLoopRange = theLoopRange;
        myLoopRangeLength = theLoopRange.length();
        N1.set(0);
    }

    /**
     * {@inheritDoc}
     *
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
     */
    public LongRange next(int theThreadIndex) {
        for (;;) {
            long oldN1 = N1.get();
            LongRange result
                    = myLoopRange.chunk(oldN1,
                            Math.max(N2, (myLoopRangeLength - oldN1) / two_K));
            long N = result.length();
            if (N == 0) {
                return null;
            }
            long newN1 = oldN1 + N;
            if (N1.compareAndSet(oldN1, newN1)) {
                return result;
            }
        }
    }

}
