//******************************************************************************
//
// File:    FixedLongSchedule.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.FixedLongSchedule
//
// This Java source file is copyright (C) 2010 by Alan Kaminsky. All rights
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

/**
 * Class FixedLongSchedule provides a fixed schedule object. The loop index is
 * type <TT>long</TT>. The loop iterations are apportioned among the parallel
 * team threads once at the beginning of the parallel for loop, with each thread
 * getting a fixed number of iterations, the same number of iterations for each
 * thread (plus or minus one).
 *
 * @author Alan Kaminsky
 * @version 27-Jan-2010
 */
class FixedLongSchedule
        extends LongSchedule {

// Hidden data members.
    // Chunk for each thread.
    private LongRange[] myChunk;

// Exported constructors.
    /**
     * Construct a new fixed schedule object.
     */
    public FixedLongSchedule() {
        super();
    }

    /**
     * Construct a new fixed schedule object. This constructor is for use by the
     * <TT>LongSchedule.parse()</TT> method.
     *
     * @param args Array of argument strings.
     */
    public FixedLongSchedule(String[] args) {
        super();
        throw new IllegalArgumentException("FixedLongSchedule(): Usage: -Dpj.schedule=fixed");
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
        return true;
    }

// Hidden operations.
    /**
     * Start generating chunks of iterations for a parallel for loop using this
     * schedule.
     * <P>
     * The <TT>start()</TT> method is only called by a single thread in the
     * Parallel Java middleware.
     *
     * @param K Number of threads in the parallel team.
     * @param theLoopRange Range of iterations for the entire parallel for loop.
     * The stride may be 1 or greater.
     */
    public void start(int K,
            LongRange theLoopRange) {
        myChunk = theLoopRange.subranges(K);
        for (int i = 0; i < K; ++i) {
            if (myChunk[i].length() == 0) {
                myChunk[i] = null;
            }
        }
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
     * @param theThreadIndex Thread index in the range 0 .. <I>K</I>-1.
     *
     * @return Chunk of iterations, or null if no more iterations.
     */
    public LongRange next(int theThreadIndex) {
        LongRange chunk = myChunk[theThreadIndex];
        myChunk[theThreadIndex] = null;
        return chunk;
    }

}
